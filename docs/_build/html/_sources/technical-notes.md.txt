# Technical Notes

## Nextflow Unstaging Hack
At the time when this pipeline was developed, Nextflow has a longstanding three-headed logic bug
in its handling of unstaging files from a scratch directory:
[nextflow-io/nextflow/issues/3711](https://github.com/nextflow-io/nextflow/issues/3711), [nextflow-io/nextflow/issues/4755](https://github.com/nextflow-io/nextflow/issues/4755), [nextflow-io/nextflow/issues/4718](https://github.com/nextflow-io/nextflow/issues/4718)

To discuss my solution we have to understand the nextflow `nxf_main` normal execution flow.

{lineno-start=1 emphasize-lines="4,6,7"}
```bash
nxf_main() {
	<...>
	set +e
	(set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
	pid=$!
	wait $pid || nxf_main_ret=$?
	nxf_unstage
	# 'afterScript' directive
}
```
the `nxf_launch` surrounded by other line noise runs the user process script, with various error
and logging options.

{lineno-start=5}
```bash
	pid=$!
	wait $pid || nxf_main_ret=$?
```
when it finishes, its exit code is saved to `nxf_main_ret`.
this variable is *the* source of the exit code.  it is used later in `on_exit()` to write to both
`.exitcode` and the nextflow wrapper's own return status.

*after* saving the exit code, `nxf_unstage` runs.  `nxf_unstage` **does not update nxf_main_ret**
if a problem occurs.  this is problem #1: nextflow does not propagate exit codes when unstaging.

problem #2 was identified by some Nextflow contributors in the issues listed above: 
in `nxf_unstage()`, `|| true` prevents stageout errors (e.g. from `nxf_fs_copy`) from propagating.
However, we actually don't mind this; it normally wouldn't be completely disastrous to continue
copying what does survive for the user to debug,
we want to survive to cleanup scratch space and run `afterScript`, et cetera.

problem #3 is that Nextflow does not provide a method to validate output after unstaging.
a 0-length `out.fastq` is just as valid as the well-formed result is just as valid as a file truncated
in the middle of the quality string.

### Manifestation
During testing, our lustre filesystem had a phase where some metadata accesses would go offsite and
fail to retrieve file information correctly.  When this (or other race conditions) occurred,
copying the output file from scratch to work directory would fail with
`cp: SRR12661676.trimmed.fastq: failed to get extents info: Invalid argument`.

In `nxf_unstage()`, each output file is being passed to `nxf_fs_copy`:

{emphasize-lines="7,9,10,11"}
```bash
nxf_unstage() {
	true
	cp .command.out	<workdir>/.command.out || true
	cp .command.err <workdir>/.command.err || true
	cp .command.trace <workdir>/.command.trace || true
	cp .command.env <workdir>/.command.env || true
	[[ ${nxf_main_ret:=0} != 0 ]] && return
	IFS=$'\n'
	for name in $(eval "ls -1d *.filtered.fastq.gz" | sort | uniq); do
		nxf_fs_copy "$name" <workdir> || true
	done
	unset IFS
}
```
{emphasize-lines="6"}
```bash
nxf_fs_copy() {
	local source=$1
	local target=$2
	local basedir=$(dirname $1)
	mkdir -p $target/$basedir
	cp -fRL $source $target/$basedir
}
```
problem 1: any error in `nxf_fs_copy` is not propagated to exitcode and hidden from the nextflow
runner.  errors are not handled at all other than (problem 2) to ignore them.

{start-line=7}
```[[ ${nxf_main_ret:=0} != 0 ]] && return```
notice the hint at problem 3: if the process script had a nonzero exit status, `nxf_unstage` avoids
unstaging whatever output files there might be, so that post-task processing doesn't see anything
that looks like valid output!

The result was that hundreds of runs would silently fail -- nextflow did not see them as errors
because one of the files in a mate pair was well-formed, or the filename existed but its contents 
could not be copied.  These runs continued through the pipeline and failed at unhelpful points,
already too late to send it back to the *true* failed step to retry.

### Solution
we need to replace nextflow's `nxf_fs_copy` with one that
 - checks for nonzero exit code from the `cp -fRL $source $target/$basedir`
 - waits before trying again, for networked filesystem race conditions
 - if still unable to unstage, *deletes outputs so that nextflow detects the failure properly*

```bash
nxf_fs_copy() {
	local source=$1
	local target=$2
	local basedir=$(dirname $1)
	mkdir -p $target/$basedir
	cp -fRL $source $target/$basedir
	nxf_fs_copy_ret=$?
	if [[ $nxf_fs_copy_ret > 0 ]]
	then
		sleep 10  # second chance
		cp -fRL $source $target/$basedir
		nxf_fs_copy_ret=$?
		if [[ $nxf_fs_copy_ret > 0 ]]
		then
			ls -alh
			echo "NXF_UNSTAGE copy $source failed ($nxf_fs_copy_ret)" | tee .command.err
			local temp=$(basename $source)
			rm -fr $target/$basedir/${temp//_[12].fastq/_*.fastq} 
		fi
	fi
}
```


## Dependencies
There are multiple ways to obtain the dependencies for the RNAquarium preprocessing pipeline
 - installing locally
   - if it works, it works.  May be unsafe for heterogeneous computing clusters, as no attempt is
     made to identify different machine types.
 - using conda packages (preferred for ease of deployment)
   - conda packages may not have optimizations for the user's hardware.
   - some tools are not available as conda packages, and `setup-minimal.sh` is provided for these.
 - using a singular "RNAquarium" container with all tools
   - may not have all optimizations for the user's hardware
   - may have limited network bandwidth, which slows down the **prefetch** step.
 - using individual biocontainer package containers
   - not all packages used are available, and `setup-minimal.sh` is provided for these.
   - may have limited network bandwidth, which slows down the **prefetch** step.

Some dependencies present additional problems which require careful testing when using a container.
newer versions of gmap/gsnap are much faster at mapping, but have a bad container, different
sensitivity, fail to detect all the input reads, or crash on certain input sequences.
```
gsnap_version,         gsnap_time,
after_dedup,total_read,aligned,multialign,aligned_unique,unaligned,mixed
2021.08.25--pl5321h...,     99.0s,
      17507,     17507,   1221,         0,          1221,    16286,    0
2021.12.17(self-compiled),  27.4s,
      17507,     17507,   1150,         0,          1150,    16357,    0
! 2023.03.24--h67092d7_0 doesn't work (bad hardcoded perl path in container)
! 2023.03.24--h9d449c0_1 doesn't work (bad hardcoded perl path in container)
! 2023.04.28--h9d449c0_1 doesn't work (bad hardcoded perl path in container)
! 2023.07.20--h9d449c0_1 doesn't work (bad hardcoded perl path in container)
2023.10.10--h9d449c0_1,     16.1s,
      17507,      2399,    158,         0,           158,     2241,    0
2023-12-01(self-compiled),   6.1s,
      17507,     17507,    275,         0,           275,    17232,    0
! 2024-02-22(self-compiled) inconsistent (segfaults on some inputs)
```
Other dependencies may or may not be safe to up- or downgrade.



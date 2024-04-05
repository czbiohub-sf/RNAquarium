include { merge_unmapped                       } from './modules/merge_unmapped.nf'
include { spades_single_end; spades_paired_end } from './modules/assembly.nf'

params.bioproj_map = "$projectDir/test/sample_mapping.json"
params.gsnap_dir   = "/hpc/projects/balla_group/sra_experiments/RNAquarium_prototyping/rnaquarium_unmapped_sample/"

bioproj_ids = Channel.fromPath(params.bioproj_map)
    .splitJson()
    .map{ it.key }

bioproj_ids = Channel.of('PRJNA834970', 'PRJNA961970', 'PRJNA828542')

workflow {
    results = merge_unmapped(bioproj_ids, params.bioproj_map, params.gsnap_dir)

    single_end_fqs = results.filter{ it[1].size() == 1 }
    paired_end_fqs = results.filter{ it[1].size() == 2 }
    both_end_fqs   = results.filter{ it[1].size() == 3 }

    // https://stackoverflow.com/a/75248731
    single_end_fqs = single_end_fqs
        .concat(
            both_end_fqs.map{ x, y -> tuple(x, y.findAll{ it =~ /PRJNA\d+_S/ }) }
        )
        .map{ id, fqs -> tuple(id, fqs[0]) }

    paired_end_fqs = paired_end_fqs
        .concat(
            both_end_fqs.map{ x, y -> tuple(x, y.findAll{ it =~ /PRJNA\d+_P/ }) }
        )
        .map{ id, fqs -> tuple(id, fqs[0], fqs[1]) }

    spades_single_end(single_end_fqs)
    spades_paired_end(paired_end_fqs)
}

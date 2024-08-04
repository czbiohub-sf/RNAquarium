import plotly.graph_objects as go
import pandas as pd

# starting_runs,successful,download_fail,corrupt_file,dedup_fail,gsnap_fail,other,preempted,dropout,unsuccessful
d = pd.read_csv("sankey.csv", header=0).to_numpy()[0]

NODES = dict( #    0                 1                          2        3       4           5
    label = [f"Starting Accessions\n({d[0]})", # 0
             f"Failed download\n({d[2]})", # 1
             f"Corrupt file?\n({d[3]})", # 2
             f"Failed GSnap\n({d[5]})", # 3
             f"Other\n({d[6]})", # 4
             f"Run 1 preemptions?\n({d[7]})", # 5
             f"Dedup dropout (too short)\n({d[4]})", # 6
             f"Successful\n({d[1]})", # 7
             f"Unsuccessful\n({d[9]})", # 8
             f"Dropout\n({d[8]-d[9]})", # 9
             ],
    pad = 15,
    thickness = 20,
    )
LINKS = dict(
    source = [ 0,  0, 0, 0, 0, 0, 8, 8, 8, 8, 8,  0, 9 ], # The origin or source nodes of
    target = [ 7,  8, 8, 8, 8, 8, 1, 2, 3, 4, 5,  9, 6 ], # The dest or target
    # d.download_fail, d.corrupt_file, d.gsnap_fail, d.other, d.preempted, d.dedup_fail,
    # d.successful, d.download_fail, d.corrupt_file, d.gsnap_fail, d.other, d.preempted,
    # d.dedup_fail
    # success,
    # download, corrupt, gsnap, other, preempted
    # download, corrupt, gsnap, other, preempted
    # 
    value = [ d[1],
              d[2], d[3], d[5], d[6], d[7],
              d[2], d[3], d[5], d[6], d[7],
              d[4], d[4] ],
    # Color of the links
    # color = ["seagreen","dodgerblue","orange", "gold", "silver","brown" ],
)
data = go.Sankey(node = NODES, link = LINKS)
fig = go.Figure(data)
fig.write_image("sankey.png")

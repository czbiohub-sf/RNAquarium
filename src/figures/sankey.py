import plotly.graph_objects as go

NODES = dict( #    0                 1                          2        3       4           5
    label = ["Starting Accessions\n(60912)", # 0
             "Successful\n(58935)", # 1
             "Failed GSnap\n(32)", # 2
             "Run 1 preemptions?\n(1371)", # 3
             "Unsuccessful\n(1945)", # 4
             "Dedup dropout (too short)", # 5
             "Failed download", # 6
             "Other", # 7
             "Final\n(60338)", # 8
             "Dropout\n(574)", # 9
             ],
    pad = 15,
    thickness = 20,
    ) # color = ["seagreen","dodgerblue","orange", "gold", "silver","brown" ],
LINKS = dict(   source = [     0,  0,    0,  4,    4,  4,  4, 1, 2, 3, 5, 6, 7 ], # The origin or the source nodes of
                target = [     1,  2,    4,  3,    5,  6,  7, 8, 8, 8, 9, 9, 9 ], # The destination or the target
                value =  [ 58935, 32, 1945, 1371, 112, 57, 405, 58935, 32, 1371, 112, 57, 405 ], # The width (quantity) of the links
                # Color of the links
             )
data = go.Sankey(node = NODES, link = LINKS)
fig = go.Figure(data)
fig.write_image("sankey.png")

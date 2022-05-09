import json
import sys

time_point = ["E3","E3.5","E4.5","E5.25","E5.5","E6.25","E6.5","E6.75","E7","E7.25","E7.5","E7.75","E8","E8.25","E8.5a","E8.5b","E9.5","E10.5","E11.5","E12.5","E13.5"]
time_point_n = len(time_point)
time_point_id = {}
for i in range(0,time_point_n):
    time_point_id[time_point[i]] = i


### read edge info 
edge = {}; node_all = []; node_each = {}; main_edge = set()

file = open("./mouse_edge.txt")
for line in file:
    l = line.rstrip().split('\t')

    if l[0] not in node_all:
        node_all.append(l[0])
    if l[1] not in node_all:
        node_all.append(l[1])

    edge[l[0]] = edge.get(l[0], [])
    edge[l[0]].append(l[1])

    num = time_point_id[l[0].split(':')[0]]
    node_each[num] = node_each.get(num, [])
    if l[0] not in node_each[num]:
        node_each[num].append(l[0])

    main_edge.add((l[0],l[1]))

### read edge weight file
main_weight = {}; extra_weight = {}
file = open("./edge_prob.txt")
for line in file:
    l = line.rstrip().split('\t')
    if float(l[2])>0.8:
        xx = "specific"
    elif float(l[2])>0.2:
        xx = "additional"
    if (l[0],l[1]) in main_edge:
        if float(l[2]) > 0.2:
            main_weight[l[1]] = float(l[2])
        else:
            main_weight[l[1]] = 0
    else:
        if float(l[2]) > 0.2:
            extra_weight[(l[0],l[1])] = float(l[2])


### read group information, used for the color of the node
file = open("./celltype_group.txt")
i = 1
coor = {}
node_group = {}
for line in file:
    l = line.rstrip().split("\t")
    node_group[l[0]] = int(l[1])
    coor[l[0]] = i
    i += 1
file.close()

### create the info for all the node

dat = {}
for i in node_all:

    dat[i] = {'name':i}
    if i in edge:
        dat[i]['children'] = []

    if i in main_weight:
        dat[i]['edge_weight'] = main_weight[i]
    else:
        print(i)

    ### 10. coors of x (used to create the plot, the cell type order)
    if i.split(':')[1] in coor:
        dat[i]['fx'] = str(coor[i.split(':')[1]])
    else:
        print(i)

    ### 11. group of celltypes. 1: neuroectoderm; 2: surface ectoderm; 3: endoderm; 4: mesoderm; 5: exe germ layer; 0: other
    ### 1: green, 2: blue, 3: red; 4: yellow; 5: brown; 0: black
    color_map = {1:"#145A32",2:"#78281F",3:"#9A7D0A",4:"#2471A3",5:"#E67E22",0:"#000000"}
    if i.split(':')[1] in node_group:
        dat[i]['node_group'] = color_map[node_group[i.split(':')[1]]]
    else:
        print(i)


### add extra node first
for i in extra_weight:
    tmp = {'name': i[1], 'edge_weight': extra_weight[i], 'fx': str(coor[i[1].split(':')[1]]), 'node_group': color_map[node_group[i[1].split(':')[1]]]}
    dat[i[0]]['children'] = dat[i[0]].get('children',[])
    dat[i[0]]['children'].append(tmp)
    print(dat[i[0]]['children'])



### connect each time point
for i in range(time_point_n-2, 0, -1):
    for j in node_each[i]:
        for k in edge[j]:
            dat[j]['children'].append(dat[k])

for k in edge['E3:Morula']:
    dat['E3:Morula']['children'].append(dat[k])
dat_json = dat['E3:Morula']


with open("./mouse.json", 'w') as json_file:
    json.dump(dat_json, json_file)



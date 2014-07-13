import sys
s = []
s.append("Metadata_position,")
#s.append(Metadata_cycle,")
for cycle in range(1,5):
    for type in range(1,5):
        s.append("Image_PathName_%d_%d,Image_FileName_%d_%d," % (cycle, type, cycle, type))

s.append("Image_PathName_DO,Image_FileName_DO\n")

for position in range(1,17):
    s.append("%d," % (position))
    for cycle in range(1,5):
        #s.append("%d," % cycle)
        # must correspond to ID_list_BCpanel.m
        for type in ["T","G","C","A"]:
            s.append("%d/%d/,%s.tif," % (position, cycle, type))
    s.append("%d/,%s.tif" % (position, "do"))
    s.append("\n")

print "".join(s)

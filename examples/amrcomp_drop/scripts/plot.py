import yt
# Load native AMReX plotfile (point to the directory)
ds = yt.load("amrviz/drop/nga2.cell.000001")
# 2D slice through z=0
p = yt.SlicePlot(ds, "z", "VF")
p.save("VF.png")
# Or multiple fields
for field in ["VF", "RHOL", "RHOG"]:
    p = yt.SlicePlot(ds, "z", field)
    p.save(f"{field}.png")

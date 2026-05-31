# ParaView Programmable Source — moving wall visualization
#
# Renders a solid cylinder (axis along x) whose right face tracks the
# analytical wall position xw = Xw + Uw*t. The wall itself is not written
# to the viz output by NGA2; this source builds it inside ParaView from
# the animation time alone.
#
# Setup in ParaView:
#   1. Sources -> Programmable Source
#   2. Properties panel: set "Output DataSet Type" to vtkPolyData
#   3. Paste the SCRIPT block below into the "Script" field
#   4. Paste the REQUEST INFORMATION SCRIPT block below into the
#      "RequestInformation Script" field
#   5. Apply. Add a Surface representation; the wall will track time.
#
# Edit Xw, Uw to match the input file's "Wall location" / "Wall velocity".
# Edit t_max to anything >= the latest time you'll scrub to (1e6 also works).


# =============================================================
# SCRIPT field
# =============================================================
import vtk

# === EDIT TO MATCH YOUR INPUT FILE ===
Xw = -3.0      # initial wall x-position
Uw = 1.0       # wall velocity
# =====================================

# Cylinder dimensions
radius = 10
length = 50.0

# Get current animation time
info = self.GetOutputInformation(0)
key = vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()
t = info.Get(key) if info.Has(key) else 0.0
xw = Xw + Uw * t

# Build a capped cylinder along default y-axis
src = vtk.vtkCylinderSource()
src.SetRadius(radius)
src.SetHeight(length)
src.SetResolution(64)
src.SetCapping(True)
src.Update()

# Rotate y-axis -> x-axis, then translate so right face is at xw
trans = vtk.vtkTransform()
trans.Translate(xw - 0.5*length, 0.0, 0.0)
trans.RotateZ(-90.0)

tf = vtk.vtkTransformPolyDataFilter()
tf.SetTransform(trans)
tf.SetInputConnection(src.GetOutputPort())
tf.Update()

output.ShallowCopy(tf.GetOutput())


# =============================================================
# REQUEST INFORMATION SCRIPT field
# =============================================================
import vtk

# === EDIT TO MATCH YOUR INPUT FILE ===
t_max = 100.0   # match Max time
# =====================================

outInfo = self.GetOutputInformation(0)
key_range = vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE()
outInfo.Set(key_range, [0.0, t_max], 2)

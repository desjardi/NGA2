# Import libraries
import os
from paraview.simple import *

desired_times = [0.000, 0.125, 0.250, 0.500]
n_times = len(desired_times)
tolerance = 1e-4
line_res = 50
res_path = './centerline'
if not os.path.exists(res_path):
    try:
        os.mkdir(res_path)
    except OSError as error:
        print(f"Error: {error}")

# Find source
source = FindSource('nga.case')

# Get the TimeKeeper object
timeKeeper = GetTimeKeeper()

# Get the number of time steps
num_time_steps = len(timeKeeper.TimestepValues)

# Get the field
T = source.CellData.GetArray('T')

# Index of the plot over line at y = 0
y0_ind = -1

# Zl
Zl_y0 = []

# Iterate over each time step
for time_step in range(num_time_steps):

    # Set the current time step
    SetActiveView(GetRenderView())
    animationScene = GetAnimationScene()
    animationScene.AnimationTime = timeKeeper.TimestepValues[time_step]

    # Get the current time
    time = animationScene.AnimationTime

    # Check if the current time is one of the specified times
    if any(abs(time - desired_time) <= tolerance for desired_time in desired_times):

        print(time)

        # Increment the index
        y0_ind = y0_ind + 1

        # create a new 'Plot Over Line'
        y0_line = PlotOverLine(registrationName=f"{'y0_'}{y0_ind}", Input=source)
        y0_line.Resolution = line_res
        y0_line.Point1 = [-0.5, 0.0, 0.0]
        y0_line.Point2 = [+0.5, 0.0, 0.0]

        # Get the VTK object associated with the plot over line
        vtk_object = y0_line.GetClientSideObject()

        # Force the update of the pipeline
        vtk_object.Update()

        # Get the VTK output data object from the VTK object
        output_data = vtk_object.GetOutput()

        # Get the CellData
        point_data = output_data.GetPointData()
        
        # Get the arrays from CellData
        Zl_y0 = point_data.GetArray('Zl')
        x_y0  = point_data.GetArray('arc_length')

        # Write to file
        with open(res_path + '/numerical_' + str(desired_times[y0_ind]) + '.dat', 'w') as file_y0:
            file_y0.write(f"x            Zl\n")
            for i in range(0, line_res):
                file_y0.write(f"{x_y0.GetValue(i)-0.5}            {Zl_y0.GetValue(i)}\n")
        file_y0.close()
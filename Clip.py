from paraview.simple import *

# Step 1: Load and merge blocks
src = GetActiveSource()
merged = MergeBlocks(Input=src)

# Step 2: Clip the left half of the cylinder (X < 0)
clip = Clip(Input=merged)
clip.ClipType = 'Plane'
clip.ClipType.Origin = [0, 0, 0]
clip.ClipType.Normal = [-1, -1, 0]
clip.Invert = 1

# Step 3: Hide original, show clip
Hide(src)
Show(clip)

# Step 4: View and display properties
view = GetActiveViewOrCreate('RenderView')
clipDisplay = GetDisplayProperties(clip, view=view)

# Step 5: Geometry info
info = clip.GetDataInformation()
print("✅ Geometry Info:")
print("  - Number of cells:", info.GetNumberOfCells())
print("  - Number of points:", info.GetNumberOfPoints())
if info.GetNumberOfCells() == 0:
    raise RuntimeError("❌ No geometry to render. Check the clip or data source.")

# Step 6: Color by K_sat_original if available
assigned_array = "K_sat_original"
pointData = info.GetPointDataInformation()
cellData = info.GetCellDataInformation()
array_names = {pointData.GetArrayInformation(i).GetName() for i in range(pointData.GetNumberOfArrays())}
array_names.update(cellData.GetArrayInformation(i).GetName() for i in range(cellData.GetNumberOfArrays()))

if assigned_array in array_names:
    ColorBy(clipDisplay, ('POINTS', assigned_array))
    clipDisplay.RescaleTransferFunctionToDataRange(True, False)
    clipDisplay.SetScalarBarVisibility(view, True)
    print(f"🎨 Coloring by: {assigned_array} (POINTS)")
else:
    print(f"⚠️ '{assigned_array}' not found. Using solid color.")
    ColorBy(clipDisplay, None)

# Step 7: Appearance
clipDisplay.Representation = 'Surface'
clipDisplay.Opacity = 1.0

# ✅ CLEAN WHITE BACKGROUND SETTINGS - FORCE WHITE
view.Background = [1.0, 1.0, 1.0]  # Pure white background
view.Background2 = [1.0, 1.0, 1.0] # Secondary background also white

# Force single color background mode
try:
    view.BackgroundColorMode = 'Single Color'
except AttributeError:
    pass

# Additional background settings
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass

# Hide ALL UI elements for cleanest output
view.OrientationAxesVisibility = 0  # Hide orientation axes
view.CenterAxesVisibility = 0       # Hide center axes

# Hide 3D widget
try: 
    clip.Show3DWidget = 0
except AttributeError: 
    pass

# ✅ DISABLE AXIS GRID for cleanest look (uncomment if you want no grid)
#view.AxesGrid.Visibility = 0

# OR if you want to keep the grid, make it black on white:
view.AxesGrid.Visibility = 1
view.AxesGrid.XTitleColor = [0, 0, 0]
view.AxesGrid.YTitleColor = [0, 0, 0]
view.AxesGrid.ZTitleColor = [0, 0, 0]
view.AxesGrid.XLabelColor = [0, 0, 0]
view.AxesGrid.YLabelColor = [0, 0, 0]
view.AxesGrid.ZLabelColor = [0, 0, 0]
view.AxesGrid.GridColor = [0.3, 0.3, 0.3]  # Light gray grid lines
view.AxesGrid.XTitleFontSize = 13          # Increased from 12
view.AxesGrid.YTitleFontSize = 13          # Increased from 12
view.AxesGrid.ZTitleFontSize = 13          # Increased from 12
view.AxesGrid.XLabelFontSize = 11          # Increased from 10
view.AxesGrid.YLabelFontSize = 11          # Increased from 10
view.AxesGrid.ZLabelFontSize = 11          # Increased from 10

# Force axis labels to be visible - ESSENTIAL
try:
    view.AxesGrid.XAxisLabels = 1
    view.AxesGrid.YAxisLabels = 1
    view.AxesGrid.ZAxisLabels = 1
    view.AxesGrid.ShowGrid = 1
    view.AxesGrid.ShowEdges = 1
    view.AxesGrid.ShowTicks = 1
    
    # Fix axis scaling - use actual data coordinates
    view.AxesGrid.XAxisUseCustomBounds = 0  # Use data bounds, not custom
    view.AxesGrid.YAxisUseCustomBounds = 0
    view.AxesGrid.ZAxisUseCustomBounds = 0
    
    # Disable axis normalization
    view.AxesGrid.DataScale = [1.0, 1.0, 1.0]  # No scaling
    view.AxesGrid.DataPosition = [0.0, 0.0, 0.0]  # No offset
    
except AttributeError:
    pass

# Make axis lines thicker
try:
    view.AxesGrid.XAxisLineWidth = 2       # Default is usually 1
    view.AxesGrid.YAxisLineWidth = 2
    view.AxesGrid.ZAxisLineWidth = 2
except AttributeError:
    # Fallback for different ParaView versions
    try:
        view.AxesGrid.LineWidth = 2
    except AttributeError:
        pass

# REMOVE ALL CUSTOM TICK SPACING - let ParaView use defaults
# Don't set any custom intervals to avoid breaking label generation

# Debug and fix axis bounds - GET ACTUAL COORDINATES
try:
    # Get the actual data bounds
    bounds = clip.GetDataInformation().GetBounds()
    print(f"Data bounds: X=[{bounds[0]:.1f}, {bounds[1]:.1f}], Y=[{bounds[2]:.1f}, {bounds[3]:.1f}], Z=[{bounds[4]:.1f}, {bounds[5]:.1f}]")
    
    # Set custom bounds using the method that works
    try:
        view.AxesGrid.CustomBounds = [bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]]
        view.AxesGrid.UseCustomBounds = 1
        print("✅ Set custom bounds using CustomBounds")
    except AttributeError:
        print("❌ Could not set custom bounds")
    
except (AttributeError, Exception) as e:
    print(f"Error getting bounds: {e}")

# Step 9: Scalar bar styling (black text on white background)
if assigned_array in array_names:
    sb = GetScalarBar(clipDisplay.LookupTable, view)
    sb.TitleColor = [0, 0, 0]      # Black title
    sb.LabelColor = [0, 0, 0]      # Black labels
    sb.TitleFontSize = 13          # Increased from 12
    sb.LabelFontSize = 11          # Increased from 10
    # Position scalar bar nicely
    sb.Position = [0.85, 0.2]      # Right side of image
    sb.ScalarBarLength = 0.6

# ✅ ADDITIONAL CLEAN OUTPUT SETTINGS
# Remove any potential gradients or lighting artifacts
try:
    view.BackgroundColorMode = 'Single Color'  # Modern ParaView versions
except AttributeError:
    try:
        view.UseGradientBackground = 0         # Older ParaView versions
    except AttributeError:
        pass

view.Background2 = [1.0, 1.0, 1.0] # Set secondary background to white too

# Improve rendering quality
try:
    view.EnableRayTracing = 0          # Disable for faster, cleaner renders
except AttributeError:
    pass

try:
    view.OSPRayMaterialLibrary = None
except AttributeError:
    pass

# Step 10: Final render and save with high quality settings
# FORCE one more background update before rendering
view.Background = [1.0, 1.0, 1.0]
Render()

# Or hide all widgets in the view
try:
    view.InteractionMode = '3D'
    view.Show3DWidgets = 0
except AttributeError:
    pass

# ✅ ENHANCED SAVE SETTINGS for cleanest PNG
try:
    SaveScreenshot("/home/arash/Projects/VN Drywell_Models/clip_clean_output.png", 
                   view,
                   ImageResolution=[1920, 1080],
                   TransparentBackground=False,  # Solid white background
                   ImageQuality=100,             # Maximum quality
                   CompressionLevel='0')         # No compression artifacts
except TypeError:
    # Fallback for older ParaView versions that don't support all parameters
    SaveScreenshot("/home/arash/Projects/VN Drywell_Models/clip_clean_output.png", 
                   view,
                   ImageResolution=[1920, 1080],
                   TransparentBackground=False)

print("✅ Clean white background screenshot saved successfully.")

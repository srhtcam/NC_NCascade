# Narrowing-Cascade Splines

![NC nets in nature](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Bunny.jpg?raw=true)

Quad-dominant meshes are popular with animation designers and can efficiently be
generated from point clouds. To join primary regions, quad-dominant meshes include
non-4-valent vertices and non-quad regions. To transition between regions of rich detail
and simple shape, quad-dominant meshes commonly use a cascade of n − 1 triangles
that reduce the number of parallel quad strips from n + 1 to 2.
For these cascades, the Narrowing-Cascade spline, short NCn, provides a new shape-
optimized G1 spline surface. NCn can treat cascade meshes as B-spline-like control
nets. For n > 3, as opposed to n = 2, 3, cascades have interior points that both guide and
complicate the construction of the output tensor-product NC spline. The NCn spline
follows the input mesh, including interior points, and delivers a high-quality curved
surface of low degree.


[Paper](https://github.com/srhtcam/NC_NCascade/raw/master/resources/25nCascade_Clean.pdf)   ,    [Cite](https://doi.org/10.1016/j.cag.2025.104441)



This code is designed to run within Grasshopper for Rhinoceros 3D to visualize surface outputs, 
meshes, and generate [BezierView](https://www.cise.ufl.edu/research/SurfLab/bview/) files (.bv) for use in this paper.

Additionally, a MATLAB script that generates NC4 surface patches for n = 4 narrowing cascades is included for users who prefer not to use Rhinoceros 3D.

## <ins>MATLAB Script for n = 4 Narrowing Cascade </ins>
This MATLAB script generates NC4 surface patches for "**n = 4**" narrowing-cascade size from a given set of control points, without the need for Rhino 3D or Grasshopper.
For the case n = 4, there are a total of 52 control points, as shown in the figure.

<p align="center">
<img src="https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/txt_pts_order.jpg?raw=true" width="300"/>
</p>

The control points should be written line by line to the .txt file, following a horizontal sequence from either the left-bottom or 
right-bottom corner upward, with no spaces between the x, y, and z coordinates, as illustrated below:

 <pre>
left-bottom to up;         right-bottom to up;
 1x,1y,1z                  8x,8y,8z
 2x,2y,2z                  7x,7y,7z
 3x,3y,3z                  6x,6y,6z
 4x,4y,4z                  5x,5y,5z
 .                         .
 .                         .
 .                         .
 51x,51y,51z			   49x,49y,49z
 52x,52y,52z			   48x,48y,48z
</pre>


When the script (*NC4_SurfPatch.m*) is executed, it prompts the user to select a *.txt* file containing the control points, and then outputs a corresponding *.bv* file in the same directory.
This .bv file can be freely visualized using [BezierView](https://www.cise.ufl.edu/research/SurfLab/bview/). 
The control points used in the paper’s figures and the MATLAB script can be downloaded from the link below:

[Download matlab.zip](https://github.com/srhtcam/NC_NCascade/raw/master/resources/matlab.zip)

In the figure below, the generated surface patch (from Fig10_b_NC4.txt) and the highlight lines can be visualized using [BezierView WebGL](https://www.cise.ufl.edu/research/SurfLab/bview/webgl/).
The resulting .bv file (Fig10_b_NC4.bv) includes both the extended net (yellow) and the NC4 surface patch (green).

<p align="center">
<img src="https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/bv_hld.png?raw=true" width="900"/>
</p>


## <ins>Installation without Compiling (Rhinoceros 3D)</ins>
To begin using Narrowing-Cascade Splines, first load the pre-compiled Grasshopper libraries, 
after ensuring Rhinoceros 3D (Rhino 8) is installed on your system. These libraries enable 
you to easily experiment with and visualize splines within the Rhino environment without needing 
to compile anything.

Inside Grasshopper, you can interactively adjust the control points of the splines. 
As you modify the points, the surface updates in real-time, allowing you to see how 
the shape evolves. Additionally, you can monitor the highlight-line distributions to evaluate 
the spline's quality and behavior.

Feel free to experiment with different control point configurations to deepen your understanding 
of the Narrowing-Cascade process and its impact on the generated surfaces.

---

### Downloading the Rhinoceros3D Sofware
This code is designed to run within Grasshopper for Rhinoceros 3D to visualize surface 
outputs, meshes, and generate Bezier view files (.bv files).

To get started, you should first install Rhinoceros 3D. You can download the software from 
the official Rhino website. If you don’t have a license, Rhino offers a 90-day free trial 
for new users. This trial allows you to fully explore and test the software's features 
before committing to a purchase.

[Download Rhino 8](https://www.rhino3d.com/download/)

---

### Pre-compiled Grasshopper Files (.gha and .gh)

You can download the pre-compiled library file (NC_NCascade.gha) and the Grasshopper file (test.gh) from the links below.

[NC_NCascade.gha](https://github.com/srhtcam/NC_NCascade/raw/master/resources/NC_NCascade.gha)

[test.gh](https://github.com/srhtcam/NC_NCascade/raw/master/NC_NCascade/test.gh)

After downloading, you need to unblock the NC_NCascade.gha file. Right-click the file, select Properties, then click Unlock, as shown below.

![Unblock the gha file](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture4.png?raw=true)

---

To install the library file in Rhino 8, first open Rhino 8 and launch Grasshopper as shown in the image below, or type the **Grasshopper** command directly into the Rhino command line.

![Launch Grasshoper](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture1.png?raw=true)

---

Inside Grasshopper, go to **File > Special Folders > Components Folder** to open the Component library file.

![Component Folder](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture2.png?raw=true)

---

Copy the downloaded NC_NCascade.gha file into the opened *component folder*. Ensure the file is checked and unblocked if needed.

![Copy .gha file](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture3.png?raw=true)


> [!TIP]
>  To install the library file (.gha), you can directly open the library folder using Windows Run.
> 1. Press *Windows + R* to open the Run dialog.
> 2. Copy and paste this path:*%Appdata%\Grasshopper\Libraries* and click *Run*.
> 3. Copy the downloaded *NC_NCascade.gha* file into this folder.

---

Restart Rhino 8 and launch Grasshopper. Then, drag the **test.gh** file onto the Grasshopper canvas to open it.

![Open the test.gh file](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture5.png?raw=true)

---

The **test.gh** file generates the default <ins>n=4</ins> Narrowing-Cascade Splines (Figure 10(b) in [paper](https://github.com/srhtcam/NC_NCascade/raw/master/resources/25nCascade_Clean.pdf)) inside Grasshopper. 
The red-colored surface in the Rhino 8 screen represents the Grasshopper surface. To use commands like Zebra Strip Analysis 
or any other Rhino 8 command, you first need to **Bake** the Grasshopper surface to convert it into a Rhino surface.

To do this, **right-click** on the **Bake this** object, select the **Bake** command, and then click **OK** in the dialog that appears.

![Bake](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture6.png?raw=true)

>Note: You can use the Bake command for all geometric objects (Geo) individually.
>The *'Bake this'* component specifically contains the outer extended surfaces and cap surfaces.

---

The **Zebra Strip Analysis** can be found under the Surface section or called directly using the command line. 
The black surface displayed over the red-colored surface represents the Rhino 8 surface generated after using the Bake command. 

![Zebra Strip Location](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture7.png?raw=true)

---

After clicking the Zebra Strip Analysis command and selecting the black surface 
(*you will notice that the red-colored surface cannot be selected*), the zebra stripes will be generated. 
To hide the red-colored surface, **right-click** on it and select **Preview**.
To show it again, right-click on the surface and select Preview once more.

![Bake](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture8.png?raw=true)

---

The **test.gh** file is briefly explained below:

![Bake](https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture9.png?raw=true)


1. Control Net Size: Adjusts the size of the control net.

2. Modifies the deformation of nets.
	-  Case: 0  -> on spherical surface (without deformation)
    -  Case: 6  -> Figure 16(c)
    -  Case: 9  -> Figure 16(e)
    -  Case: 21 -> Figure 16(g)
    -  Case: 26 -> Figure 18(a)
    -  Case: 3  -> Figure 18(c)
    -  Case: 25 -> Figure 18(d)
    -  Case: 22 -> Figure 18(e)
    -  Case: 13 -> Figure 18(f)

3. Displays mesh examples from Figure 2 (in the paper) along with their corresponding Narrowing-Cascade nets.
To view the results, connect the *..._pts* files to the *NC4_NCascade* object.

<p align="center">
<img src="https://github.com/srhtcam/NC_NCascade/blob/master/resources/images/Picture10.png?raw=true" width="900"/>
</p> 

4. Clicking the button generates the .bv file, which is saved in the *C:\Users\ ...\Downloads* directory and named surface.bv for the current surface.



## <ins>Installation with Compiling (Rhinoceros 3D) </ins>

Please install the necessary tools for Rhinoceros Scripting on Windows by following the tutorial steps below.

[Installing Tools (Windows)](https://developer.rhino3d.com/guides/rhinocommon/installing-tools-windows/)

After that, clone this repository using Visual Studio and build the project. 
This project is built for Rhino 7, so Rhino 8 users may need to adjust the dependencies, 
such as the Grasshopper version and RhinoCommon version, based on the Rhino 8 version.


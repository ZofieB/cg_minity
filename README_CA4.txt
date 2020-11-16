This implementation of the first project includes an scene setup through ray-object-intersections and a raytracing attempt
---Controls---
In the program, the user can use the inserted controls to modify the used values.
Reflection constant:
	used to control the impact of reflection in the raytracing model
Illumination constant:
	used to control the impact of the illumination model in the raytracing model
Refraction constant:
	used to control the impact of refraction in the raytracing model
Enable RT-scene:
	used to make the scene-setup of the raytracing assignment visible
Model Enabled (in menu 'Model'):
	used to enable/disable the loaded model in the viewer

---Functionality---
Intersections:
	The assignment supports intersection tests for boxes, spheres, cylinders and planes.
	The intersection tests for boxes and cylinders also support transformation of the primitives defined by matrices.
Scene setup:
	The scene can be constructed as desired using the given primitive structures. For adding or removing primitives, one has to
	modify the scene struct and add/remove the desired primitives to it.
Raytracing:
	The implemented raytracing is based on the idea of limited amount of bounces (here 4).
	A ray is first traced fully, saving all the bounced rays in an array and after all rays are traced, the final color value is 
	calculated out of all the bounced rays.

For the raytracing part as well as the solution of quadric equations for torus intersection I partly used code I found online, that I modified for my needs. 
I couldn't make torus intersection work, but I kept the code in to show what I tried. From my general understanding the intersection is 
the same as with other primitives and the hard point is in calculating the solution to the quadric equation.
Raytracing doesn't work perfectly. As I had huge issues with performance, I didn't get to debug the program using the visuals (even with a 
downscaled window). There are some reflection in the scene, but also huge noise artifacts. I suspect they are due to some floating-point-accuracy-issue, which I also 
encountered with the calculation of box-normals, but I couldn't find the errors. I generally had big issues with this part, although I would say I understand raytracing itself pretty well.
I look forward on the feedback on this :)
To work on the scene itself I commented the raytracing part out and used only the illumination model without reflection or refraction. This is also stated in the code.

This implementation of the first project includes an exploded view and animation with interpolation.
---Controls---
In the program, the user can use the inserted controls to modify the used values.
explosion factor slider
	the user can control the exploded view with this slider
add keyframe: a
	with pressing a, the user can add a keyframe for the animation
remove keyframe: r
	with pressing r, the user can remove te latest keyframe from the animation
play animation: p
	with pressing p, the user can play the animation from the 4 saved keyframes

---Functionality---
Exploded view:
	The exploded view is implemented by first computing the center of each vertex group. This is done while 
	loading the model. Each group center is stored in the group struct and therefore accessible in the ModelRenderer.cpp
	The group center for each vertex as well as the explosion factor retrieved from the GUI are passed to 
	the vertexshader for the computation of the new position of each vertex
Cubic Bezier interpolation:
	The interpolation is done for scalars, vectors and matrices seperately. Each function is located in the 
	CameraInteractor.cpp as a member function of this class. Scalar and vector interpolation are a direct implementation 
	of the cubic bezier interpolation formula.
	Matrix interpolation is achieved by decomposing the matrix first, with the use of some help functions. Each 
	part is then interpolated seperately as a vector (or quaternion) and transformed into a matrix again. The 
	three resulting matrices are then multiplied together to form the interpolated matrix.
	Translation and sclaing are simply interpolated with vector interpolation and the quaternions are SLERP interpolated 
	using the glm function mix().
Animation:
	Adding and removing a keyframe as well as playing the animation are initialized in the CameraInteractor as well. 
	Pressing a key starts the corresponding routine.
	Adding and removing a keyframe corresponds to maintaining the vector containing the keyframes for the scene. 
	This keyframes vector is part of the CameraInteractor class.
	The animation routine takes the keyframe vector and for a given time interval (here 4 seconds) interpolates 
	between these keyframes each loop iteration with the current index. The index is retrieved by the current time and 
	therefore independet of the framerate.

I tried to record a video, but my laptop is very slow and minity itself is sometimes lagging, so the videos didn't 
turn out very good. The functionality is still there and you could just try it out yourself :)
This implementation of the first project includes the Phong illumination model.
In the program, the user can use the inserted controls to modify the used values.

	shininess slider:
		this slider lets the user control the used shininess as an integer in the value range of 5 to 100
		the used value is passed to the fragment shader and used there to compute the illumination model
	ambient, specular and diffuse color:
		the color pickers let the user pick the desired color for each of the given compontents of the light source
		these are passed to the fragment shader as the respective values of the light source

The following values are passed in addition to the fragment shader:
	ambient, specular and diffuse material properties:
		these properties are passed to the shader, if they are present. If they are not present, they each get assigned a default value of (1.0, 1.0, 1.0)
	normal Matrix:
		is used to transform the normals according to the transformation model before using them for the computation

With all these passed variables, the Phong illumination model is calculated in the fragment shader.
Every used vector is normalized in the final calculation. The reflection is only considered, if it is actual pointing in the viewers direction, else it is 0.
The view and light vector are computet by subtracting the fragment position from the corresponding points (light source position and view origin position).

There are no additional features implemented.
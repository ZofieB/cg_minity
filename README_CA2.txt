This implementation of the first project includes texture maps, normal mapping and bump mapping.
---Controls---
In the program, the user can use the inserted controls to modify the used values.

	texture checkboxes:
		for each different texture there is a checkbox to enable or disable the specific texture
	normal mapping menu:
		the user can pick the kind of normal mapping that should be used or disable it with choosing None.
	bump mapping checkbox:
		only if normal mapping is used, the checkbox to enable bump mapping is available
	bump mapping parameters:
		if bump mapping is enabled, the user is able to control the different parameters that change the bump map
		and choose, which bump map function should be used. There are two different periodic functions available

---Passed variables---
The following values are passed in addition to the fragment shader:
	ambient, diffuse and specular texture:
		these properties are passed to the shader, if they are present. If they are not present, the corresponding flags are false and 
		no texture is applied, even when the checkbox is set to true
	normal mapping textures:
		these properties are also passed, if they are present and controlled by two boolean flags like the textures above.
	bump mapping parameters:
		for bump mapping there are three flags passed to enable bump mapping and the two different functions. The parameters 
		are also passed from the GUI

---Functionality---
Exercise 1
	Every texture is simply sampled for the correspoding value in the fragment shader, using the texture coordinates 
	of each fragment. The value from the texture is then multiplied on the result, to add it on top.
Exercise 2
	object space normal mapping:
		the normal is retrieved by simply sampling the corresponding texture for the color value. To get the correct normal,
		the values need to be mapped to the range -1 to 1, by multiplying by 2 and subtracting 1. This normal is then used for 
		all the lighting computations instead of the passed fragment normal

	tangent space normal mapping:
		to retrieve the correct normal, the tangent space base needs to be computet beforehand in the geometry shader. Using 
		te texture coordinates in the texture space, the tangent and bitangent are computed for each vertex and combined with 
		the vertex normal they yield the TBN matri, which is then passed to the fragment shader. There the normals that are 
		sampled from the texture, need to be transformed to the correct space (from tangent space to world space) by using the 
		TBN matrix.
		>> Although I followed the theoretical instructions to tangent space normal mapping closely (and I would say that I do 
		understand the theory behind it), there seems to be something wrong with the transformation. In the result there are some 
		primitives slightly visible. I played around a lot to find the error, but the only thing I found out is, that when I 
		do not transform the normal with the TBN matrix (what should be very neccessary), the result looks perfectly good. Because 
		of that I assume that I made some mistake in loading the correct texture, as the result should be off when not transforming 
		the normal with the TBN matrix
Exercise 3
	for bump mapping I used two different sinus functions, that the user can choose. With the already computed tangent and bitangent,
	I compute the new normal according to the equations from the lecture, using the partial derivatives of the bump function.
	I am sure, that this result is not correct yet, but I also don't really know where I went wrong. I tried a lot of different things, 
	including using different normals, with and without interpolation of tangent and bitangent and other equations I found online.
	I don't feel super confident in the theory behind bump mapping, as we didn't really have it in the lecture, so I am looking forward 
	to seeing the mistakes I made :)

There are no additional features implemendet (as I already struggled with the mandatory part).
I changed the source of the shininess for the lighting, so that now, if shininess i present in the material file, the shininess from 
that file will be used instead of the shininess retrieved by the slider.
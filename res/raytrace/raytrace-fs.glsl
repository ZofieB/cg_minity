#version 400
#extension GL_ARB_shading_language_include : require
#include "/raytrace-globals.glsl"

uniform mat4 modelViewProjectionMatrix;
uniform mat4 inverseModelViewProjectionMatrix;

uniform vec3 worldLightPosition;

uniform float currentTime;

in vec2 fragPosition;
out vec4 fragColor;

//intersection tests

struct Scene
{
	vec3 something;
};
//sphere
struct Sphere
{
	bool hit;
	vec3 near;
	vec3 far;
	vec3 normal;
};

struct Box
{
	bool hit;
	vec3 near;
	vec3 far;
};

struct Plane
{
	bool hit;
	vec3 pos; //hit point
	vec3 normal;
};

struct Cylinder
{
	bool hit;
	vec3 near;
	vec3 far;
	vec3 normal_near;
	vec3 normal_far;
};

Sphere calcSphereIntersection(float r, vec3 rayOrigin, vec3 center, vec3 rayDirection)
{
	vec3 oc = rayOrigin - center;
	vec3 dir = normalize(rayDirection);
	float dir_oc = dot(dir, oc);
	float under_square_root = dir_oc * dir_oc - dot(oc, oc) + r * r;
	if(under_square_root > 0.0)
	{
		float da = -dir_oc + sqrt(under_square_root);
		float ds = -dir_oc - sqrt(under_square_root);
		vec3 near = rayOrigin +  min(da, ds) * dir;
		vec3 far = rayOrigin + max(da, ds) * dir;
		vec3 normal = (near - center);

		return Sphere(true, near, far, normal);
	}
	else
	{
		return Sphere(false, vec3(0), vec3(0), vec3(0));
	}

}
//box -> non-oriented
//oriented: rotate box OR rotate ray -> include transformtion matrix and transform raydirection and rayorigin by inverse transformation
Box calcBoxIntersection(vec3 rayOrigin, vec3 rayDirection, vec3 boxMin, vec3 boxMax, mat4 transformation)
{
	rayOrigin = vec3(inverse(transformation) * vec4(rayOrigin, 1.0));
	rayDirection = vec3(inverse(transformation) * vec4(rayDirection, 1.0));
	vec3 tMin = (boxMin - rayOrigin) / rayDirection;
	vec3 tMax = (boxMax - rayOrigin) / rayDirection;
	vec3 t1 = min(tMin, tMax);
	vec3 t2 = max(tMin, tMax);
	float tNear = max(max(t1.x, t1.y), t1.z);
	float tFar = min(min(t2.x, t2.y), t2.z);
	if(tNear > tFar)
	{
		return Box(false, vec3(0), vec3(0));
	}
	else
	{
		return Box(true, vec3(transformation * vec4(rayOrigin + tNear * rayDirection, 1.0)), vec3(transformation * vec4(rayOrigin + tFar * rayOrigin, 1.0)));
	}
	//to get 3D position: rayOrigin + {tNear, tFar} * rayDirection
}
//plane equation : n * (q - p) = 0
Plane calcPlaneIntersection(vec3 rayOrigin, vec3 rayDirection, vec3 n, vec3 p)
{
	vec3 normal = normalize(n);
	float counter = dot(n, (p - rayOrigin));
	float numerator = dot(n, rayDirection);
	if (numerator == 0)
	{
		return Plane(false, vec3(0.0), vec3(0.0));
	}
	float s =  counter / numerator;
	if (s < 0)
	{
		return Plane(false, vec3(0.0), vec3(0.0));
	}
	else
	{
		return Plane(true, rayOrigin + s * rayDirection, normal);
	}
}
//cylinder
Cylinder calcCylinderIntersection(vec3 rayOrigin, vec3 rayDirection, float radius, float height)
{
	//y-axis aligned cylinder
	float a = rayDirection.x * rayDirection.x + rayDirection.z * rayDirection.z;
	float b = 2 * rayOrigin.x * rayDirection.x + 2 * rayOrigin.z * rayDirection.z;
	float c = rayOrigin.x * rayOrigin.x + rayOrigin.z * rayOrigin.z - radius * radius;

	float under_square_root = b*b - 4*a*c;
	if (under_square_root > 0.0)
	{
	//calc intersection with infinite cylinder
		float t1 = -b + sqrt(b*b - 4*a*c);
		t1 = t1 / (2 * a);
		float t2 = -b - sqrt(b*b - 4*a*c);
		t2 = t2 / (2 * a);
		Plane p1 = calcPlaneIntersection(rayOrigin, rayDirection, vec3(0.0, 1.0, 0.0), vec3(0.0, 0.0, 0.0));
		Plane p2 = calcPlaneIntersection(rayOrigin, rayDirection, vec3(0.0, 1.0, 0.0), vec3(0.0, height, 0.0));
		if (max(t1, t2) < 0)
		{
			//both scalars are negative -> no intersection
			return Cylinder(false, vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f));
		}
		vec3 near = rayOrigin + min(t1, t2) * rayDirection;
		vec3 far = rayOrigin + max(t1, t2) * rayDirection;
		vec3 normal_near = normalize(near - vec3(0.0, near.y, 0.0));
		vec3 normal_far = normalize(far - vec3(0.0, far.y, 0.0));
		if(near.y > height && height > far.y)
		{
		//intersection with upper boundary
			
			float t3 = (height - rayOrigin.y) / rayDirection.y;
			return Cylinder(true, rayOrigin + t3 * rayDirection, near, vec3(0.0, 1.0, 0.0), normal_near);
		}
		if (near.y < 0 && 0 < far.y)
		{
		//intersection with lower boundary
			float t3 = (0 - rayOrigin.y) / rayDirection.y;
			return Cylinder(true, rayOrigin + t3 * rayDirection, near, vec3(0.0, 1.0, 0.0), normal_near);
		}
		//calc upper and lower boundary
		if (near.y < 0 || near.y > height)
		{
			return Cylinder(false, near, far, normal_near, normal_far);
		}
		return Cylinder(true, near, far, normal_near, normal_far);
	}
	else
	{
		return Cylinder(false, vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f));
	}
	//check for boundary planes
}


float calcDepth(vec3 pos) //pass world space position -> returns depth value for OpenGL (pass to gl_FragDepth)
{
	float far = gl_DepthRange.far; 
	float near = gl_DepthRange.near;
	vec4 clip_space_pos = modelViewProjectionMatrix * vec4(pos, 1.0);
	float ndc_depth = clip_space_pos.z / clip_space_pos.w;
	return (((far - near) * ndc_depth) + near + far) / 2.0;
}

float phongIllumination(vec3 position, vec3 rayOrigin, vec3 normal, float ka, float kd, float ks, float shininess)
{
	//direction vector from intersection point to light source
	vec3 light = normalize( worldLightPosition - position );

	//rayDirection vector from intersection point to ray origin
	vec3 viewer = normalize( rayOrigin - position );

	//calculate halfway vector
	vec3 halfway = vec3(0) ; //TODO!
	float specular = max (dot( halfway, normal), 0.0);

	//calculate illumination
	return vec4( (ka * vec3(0.5, 0.5, 0.5) + kd * dot(light, normal) + ks * pow( specular, shininess) ), 1.0f);
}
vec3 closestIntersection(Scene scene)
{
	return scene.something;
}

void main()
{
	vec4 near = inverseModelViewProjectionMatrix*vec4(fragPosition,-1.0,1.0); //world space ray position (from screen space back to world space)
	near /= near.w;	//perspective divide -> diverge rays

	vec4 far = inverseModelViewProjectionMatrix*vec4(fragPosition,1.0,1.0);
	far /= far.w;

	// this is the setup for our viewing ray
	vec3 rayOrigin = near.xyz;
	vec3 rayDirection = normalize((far-near).xyz);

	Scene scene;
	int max_bounces = 4;

	/*
	for(int i = 0; i < max_bounces; i++)
	{
		//do raytracing
		//intersection tests
		vec3 p = closestIntersection(scene);
		vec3 n;

		//illumintation - shadow ray
		fragColor += phongIllumination(p, rayOrigin, normal, ka, kd, ks, shininess);

		//reflection ray
		vec3 r = (2 * dot(n, rayOrigin - p)) * n - (rayOrigin - p);

		//refraction ray
		shadowRay = 
		//set new ray Origin and direction
		rayOrigin = ;
		rayDirection = ;
	}
	*/

	vec3 spherePos = vec3(0.0f);
	vec3 boxPos = vec3(0.0f);
	vec3 planePos = vec3(0.0f);
	vec3 planeNormal = vec3(0.0, 1.0, 0.0);
	//radius * cos/sin for circular path
	spherePos.x = 1.0 * cos(currentTime);
	spherePos.z = 1.0 * sin(currentTime);

	boxPos.x = 1.0 * cos(currentTime);
	boxPos.y = 1.0 * sin(currentTime);

	mat4 boxTransformation = mat4(1.0f);
	boxTransformation[0][0] = cos(radians(45));
	boxTransformation[1][1] = cos(radians(45));
	boxTransformation[1][0] = - sin(radians(45));
	boxTransformation[0][1] = sin(radians(45));

	Box b = calcBoxIntersection(rayOrigin, rayDirection, boxPos, boxPos + vec3(0.5, 0.5, 0.5), boxTransformation);

	Sphere s = calcSphereIntersection(0.5, rayOrigin, spherePos, rayDirection);

	Plane p = calcPlaneIntersection(rayOrigin, rayDirection, planeNormal , planePos);

	Cylinder c = calcCylinderIntersection(rayOrigin, rayDirection, 0.25, 0.5);
	fragColor = vec4(vec3(0.0f), 1.0f);
	float boxDepth = -1.0f;
	if(c.hit)
	{
		fragColor+= vec4(0.0f, 1.0f, 1.0f, 1.0f);
		gl_FragDepth = calcDepth(c.near);
		return;
	}
	/*
	if(p.hit)
	{
		fragColor+= vec4(0.0f, 1.0f, 1.0f, 1.0f);
		gl_FragDepth = calcDepth(p.pos);
		return;
	}
	if(b.hit)
	{
		fragColor += vec4(0.0f, 1.0f, 0.0f, 1.0f);
		boxDepth = calcDepth(b.near);
		gl_FragDepth = boxDepth;
		return;
	}*/

	/*if(s.hit)
	{
		float intensity = dot(-rayDirection, s.normal); //simple model: assume light at camera position
		fragColor = vec4(intensity, 0.0, 0.0, 1.0);
		float sphereDepth = calcDepth(s.near);
		gl_FragDepth = max(boxDepth, sphereDepth);
	}*/
	//fragColor = vec4(1.0);

	/*where to go now:
		1. scene representation: encapsulate scene in some kind of function -> can be static setup or user interface
			-> go through primitives to find e.g. first intersection point
		2. recursive raytracer call (no recursive calls n glsl)
			-> light interactions are additive: translate recursion to iteration and count number of bounces
			-> trace a ray and branch for different kind of rays in function
		3. creativity: distributed raytracing, anti-aliasing
			shoot more rays per fragment 
				-> for loop in main and offset on rayPosition, random position in x and y and divide by viewport size
				-> average over rays or weight according to gaussian values
	*/

	// using calcDepth, you can convert a ray position to an OpenGL z-value, so that intersections/occlusions with the
	// model geometry are handled correctly, e.g.: gl_FragDepth = calcDepth(nearestHit);
	// in case there is no intersection, you should get gl_FragDepth to 1.0, i.e., the output of the shader will be ignored

	gl_FragDepth = 1.0; //-> visibility and intersection

	/*questions:
	for basic illumination in raytracing: light source attributes ia, id, is neccessary?
	no need to take models into account -> what values for ka kd and ks?
	*/
}
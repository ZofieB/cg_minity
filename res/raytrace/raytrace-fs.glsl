#version 400
#extension GL_ARB_shading_language_include : require
#include "/raytrace-globals.glsl"

uniform mat4 modelViewProjectionMatrix;
uniform mat4 inverseModelViewProjectionMatrix;

uniform vec3 worldLightPosition;

uniform float currentTime;
uniform float maxFloat;

in vec2 fragPosition;
out vec4 fragColor;

//spline function for trajectories
vec3 cubicBezierVector(vec3 cp0, vec3 cp1, vec3 cp2, vec3 cp3, float u)
{
	return pow((1 - u), 3) * cp0 + 3 * u * pow((1 - u), 2) * cp1 + 3 * pow(u, 2) * (1 - u) * cp2 + pow(u, 3) * cp3;
}
//Structs to define the scene
struct Box 
{
	vec3 bmin;
	vec3 bmax;
	mat4 transformation;
	vec3 ka;
	vec3 kd;
	vec3 ks;
	int shininess;
};
struct Sphere 
{
	vec3 center;
	float radius;
	vec3 ka;
	vec3 kd;
	vec3 ks;
	int shininess;
};
struct Plane 
{
	vec3 normal;
	vec3 point;
	vec3 ka;
	vec3 kd;
	vec3 ks;
	int shininess;
};
struct Cylinder 
{
	float radius;
	float height;
	mat4 transformation;
	vec3 ka;
	vec3 kd;
	vec3 ks;
	int shininess;
};
struct Scene
{
	Box box;
	Sphere sphere;
	//Plane plane;
	Cylinder cylinder;
};

//Intersection structs
struct Intersection
{
	bool hit;
	float tMin;
	vec3 normal;
	vec4 color;
	vec3 ka;
	vec3 kd;
	vec3 ks;
	int shininess;
};
struct SphereIntersection
{
	bool hit;
	vec3 near;
	vec3 far;
	vec3 normal;
	float t_near;
};

struct BoxIntersection
{
	bool hit;
	vec3 near;
	vec3 far;
	vec3 normal_near;
	vec3 normal_far;
	float t_near;
};

struct PlaneIntersection
{
	bool hit;
	vec3 pos; //hit point
	float t_near;
};

struct CylinderIntersection
{
	bool hit;
	vec3 near;
	vec3 far;
	vec3 normal_near;
	vec3 normal_far;
	float t_near;
};

SphereIntersection calcSphereIntersection(float r, vec3 rayOrigin, vec3 center, vec3 rayDirection)
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

		return SphereIntersection(true, near, far, normal, min(da,ds));
	}
	else
	{
		return SphereIntersection(false, vec3(0), vec3(0), vec3(0), 0);
	}

}
//box -> non-oriented
//oriented: rotate box OR rotate ray -> include transformtion matrix and transform raydirection and rayorigin by inverse transformation
BoxIntersection calcBoxIntersection(vec3 rayOrigin, vec3 rayDirection, vec3 boxMin, vec3 boxMax, mat4 transformation)
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
		return BoxIntersection(false, vec3(0), vec3(0), vec3(0), vec3(0), 0);
	}
	else
	{
		vec3 center =  (boxMin + boxMax) * 0.5;
		vec3 tmp =  (rayOrigin + tNear * rayDirection) - center;
		vec3 d = (boxMin - boxMax) * 0.5;
		float bias = 1.0001;
		vec3 normal_near = normalize( vec3(int((tmp.x / abs(d.x)) * bias), int((tmp.y / abs(d.y)) * bias), int((tmp.z / abs(d.z)) * bias)) );
		tmp = (rayOrigin + tFar * rayOrigin) - center;
		vec3 normal_far = normalize( vec3(int(tmp.x), int(tmp.y), int(tmp.z)) );
		return BoxIntersection(true, vec3(transformation * vec4(rayOrigin + tNear * rayDirection, 1.0)), vec3(transformation * vec4(rayOrigin + tFar * rayOrigin, 1.0)), normal_near, normal_far, tNear);
	}
	//to get 3D position: rayOrigin + {tNear, tFar} * rayDirection
}
//plane equation : n * (q - p) = 0
PlaneIntersection calcPlaneIntersection(vec3 rayOrigin, vec3 rayDirection, vec3 n, vec3 p)
{
	vec3 normal = normalize(n);
	float counter = dot(n, (p - rayOrigin));
	float numerator = dot(n, rayDirection);
	if (numerator == 0)
	{
		return PlaneIntersection(false, vec3(0.0), 0);
	}
	float s =  counter / numerator;
	if (s < 0)
	{
		return PlaneIntersection(false, vec3(0.0), 0);
	}
	else
	{
		return PlaneIntersection(true, rayOrigin + s * rayDirection, s);
	}
}
//cylinder
CylinderIntersection calcCylinderIntersection(vec3 rayOrigin, vec3 rayDirection, float radius, float height)
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
		PlaneIntersection p1 = calcPlaneIntersection(rayOrigin, rayDirection, vec3(0.0, 1.0, 0.0), vec3(0.0, 0.0, 0.0));
		PlaneIntersection p2 = calcPlaneIntersection(rayOrigin, rayDirection, vec3(0.0, 1.0, 0.0), vec3(0.0, height, 0.0));
		if (max(t1, t2) < 0)
		{
			//both scalars are negative -> no intersection
			return CylinderIntersection(false, vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f), 0);
		}
		vec3 near = rayOrigin + min(t1, t2) * rayDirection;
		vec3 far = rayOrigin + max(t1, t2) * rayDirection;
		vec3 normal_near = normalize(near - vec3(0.0, near.y, 0.0));
		vec3 normal_far = normalize(far - vec3(0.0, far.y, 0.0));
		if(near.y > height && height > far.y)
		{
		//intersection with upper boundary
			
			float t3 = (height - rayOrigin.y) / rayDirection.y;
			return CylinderIntersection(true, rayOrigin + t3 * rayDirection, near, vec3(0.0, 1.0, 0.0), normal_near, t3);
		}
		if (near.y < 0 && 0 < far.y)
		{
		//intersection with lower boundary
			float t3 = (0 - rayOrigin.y) / rayDirection.y;
			return CylinderIntersection(true, rayOrigin + t3 * rayDirection, near, vec3(0.0, 1.0, 0.0), normal_near, t3);
		}
		//calc upper and lower boundary
		if (near.y < 0 || near.y > height)
		{
			return CylinderIntersection(false, near, far, normal_near, normal_far, 0);
		}
		return CylinderIntersection(true, near, far, normal_near, normal_far, min(t1,t2));
	}
	else
	{
		return CylinderIntersection(false, vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f), 0);
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

float phongIllumination(vec3 position, vec3 rayOrigin, vec3 normal, vec3 ka, vec3 kd, vec3 ks, int shininess)
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
Intersection closestIntersection(Scene scene, vec3 rayOrigin, vec3 rayDirection)
{
	//calculate intersections
	BoxIntersection b = calcBoxIntersection(rayOrigin, rayDirection, scene.box.bmin, scene.box.bmax, scene.box.transformation);
	SphereIntersection s = calcSphereIntersection(scene.sphere.radius, rayOrigin, scene.sphere.center, rayDirection);
	CylinderIntersection c = calcCylinderIntersection(rayOrigin, rayDirection, scene.cylinder.radius, scene.cylinder.height);
	//find closest intersection -> minimal t_near value, if there is one
	float tMin = maxFloat;
	int shininess;
	vec3 normal, ka, kd, ks;
	vec4 color;
	if(! (b.hit || s.hit || c.hit))
	{
		//no intersection
		return Intersection(false, tMin, normal, color, vec3(0), vec3(0), vec3(0), 0);
	}
	if (b.hit && b.t_near < tMin)
	{
		tMin = b.t_near;
		normal = b.normal_near;
		color = vec4(1.0f, 0.0f, 1.0f, 1.0f);
		ka = scene.box.ka;
		kd = scene.box.kd;
		ks = scene.box.ks;
		shininess = scene.box.shininess;
	}
	if (s.hit && s.t_near < tMin)
	{
		tMin = s.t_near;
		normal = s.normal;
		color = vec4(1.0f, 1.0f, 0.0f, 1.0f);
		ka = scene.sphere.ka;
		kd = scene.sphere.kd;
		ks = scene.sphere.ks;
		shininess = scene.box.shininess;
	}
	if (c.hit && c.t_near < tMin)
	{
		tMin = c.t_near;
		normal = c.normal_near;
		color = vec4(0.0f, 1.0f, 1.0f, 1.0f);
		ka = scene.cylinder.ka;
		kd = scene.cylinder.kd;
		ks = scene.cylinder.ks;
		shininess = scene.cylinder.shininess;
	}
	return Intersection(true, tMin, normal, color, ka, kd, ks, shininess);
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

	float idx = (sin(0.5 * currentTime) + 1) / 2.0f;
	vec3 spherePos = cubicBezierVector(vec3(0.0f), vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 1.0f, 0.0f), vec3 (0.0f), idx);
	vec3 boxPos = vec3(0.0f);
	vec3 cylinderPos = cubicBezierVector(vec3(-1.0f), vec3(-1.0f, 0.0f, 0.0f), vec3(-1.0f, -1.0f, -1.0f), vec3 (-1.0f), idx);

	mat4 boxTransformation = mat4(1.0f);
	boxTransformation[0][0] = cos(radians(45));
	boxTransformation[1][1] = cos(radians(45));
	boxTransformation[1][0] = - sin(radians(45));
	boxTransformation[0][1] = sin(radians(45));

	vec3 ka_mtl = vec3(0.0f);
	vec3 kd_mtl = vec3(0.7f);
	vec3 ks_mtl = vec3(0.2f);
	int shininess = 50;

	Box box = Box(cylinderPos, cylinderPos + 0.5, boxTransformation, ka_mtl, kd_mtl, ks_mtl, shininess);
	Sphere sphere = Sphere(spherePos, 0.5, ka_mtl, kd_mtl, ks_mtl, shininess);
	Plane plane = Plane(vec3(0.0f, 1.0f, 0.0f), vec3(-1.0f, -1.0f, -1.0f), ka_mtl, kd_mtl, ks_mtl, shininess);
	Cylinder cylinder = Cylinder(0.4, 0.7, boxTransformation, ka_mtl, kd_mtl, ks_mtl, shininess);

	Scene scene = Scene(box, sphere, cylinder);
	Intersection sec = closestIntersection(scene, rayOrigin, rayDirection);

	if(sec.hit)
	{
		fragColor = sec.color * phongIllumination(rayOrigin + sec.tMin * rayDirection, rayOrigin, sec.normal, sec.ka, sec.kd, sec.ks, sec.shininess);
		gl_FragDepth = calcDepth(rayOrigin + sec.tMin * rayDirection);
		return;
	}
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
	fragColor = vec4(1.0);
	gl_FragDepth = 1.0; //-> visibility and intersection
}
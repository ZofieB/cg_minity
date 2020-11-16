#version 400
#extension GL_ARB_shading_language_include : require
#include "/raytrace-globals.glsl"

#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif

#define     EQN_EPS     1e-9
#define	    IsZero(x)	((x) > -EQN_EPS && (x) < EQN_EPS)


uniform mat4 modelViewProjectionMatrix;
uniform mat4 inverseModelViewProjectionMatrix;

uniform vec3 worldLightPosition;
uniform bool enableRT;

uniform float currentTime;
uniform float maxFloat;
uniform float kr;
uniform float kt;
uniform float kl;
uniform float n0 = 1.0f;

in vec2 fragPosition;
out vec4 fragColor;

//spline function for trajectories
vec3 cubicBezierVector(vec3 cp0, vec3 cp1, vec3 cp2, vec3 cp3, float u)
{
	return pow((1 - u), 3) * cp0 + 3 * u * pow((1 - u), 2) * cp1 + 3 * pow(u, 2) * (1 - u) * cp2 + pow(u, 3) * cp3;
}

double cbrt(double x)
{
	float third = 1.0f/3.0f;
	if (x > 0.0f)
	{
		return pow(float(x), third);
	}
	else if (x < 0.0f)
	{
		return -pow(float(-x), third);
	}
	else
	{
		return 0.0f;
	}
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
struct Torus
{
	float major_radius;
	float minor_radius;
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
	Plane plane;
	Cylinder cylinder;
	//Torus torus;
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

struct TorusIntersection
{
	bool hit;
	vec3 near;
	vec3 normal_near;
	float t_near;
};

//ray struct
struct ray
{
	vec3 p, d, n;
	vec4 color;
	float rfl, rfr;
	float n0, n1, l;
	int level, i0, i1;
	Intersection sec;
};

int SolveQuadric(double c[3], double s[2])
{
    double p, q, D;

    /* normal form: x^2 + px + q = 0 */

    p = c[ 1 ] / (2 * c[ 2 ]);
    q = c[ 0 ] / c[ 2 ];

    D = p * p - q;

    if (IsZero(D))
    {
	s[0] = - p;
	return 1;
    }
    else if (D < 0)
    {
	return 0;
    }
    else /* if (D > 0) */
    {
	double sqrt_D = sqrt(D);

	s[0] =   sqrt_D - p;
	s[1] = - sqrt_D - p;
	return 2;
    }
}


int SolveCubic(double c[4], double s[3])
{
    int i, num;
    double sub;
    double A, B, C;
    double sq_A, p, q;
    double cb_p, D;

    /* normal form: x^3 + Ax^2 + Bx + C = 0 */

    A = c[2] / c[3];
    B = c[1] / c[3];
    C = c[0] / c[3];

    /*  substitute x = y - A/3 to eliminate quadric term:
	x^3 +px + q = 0 */

    sq_A = A * A;
    p = 1.0/3 * (- 1.0/3 * sq_A + B);
    q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);

    /* use Cardano's formula */

    cb_p = p * p * p;
    D = q * q + cb_p;

    if (IsZero(D))
    {
	if (IsZero(q)) /* one triple solution */
	{
	    s[0] = 0;
	    num = 1;
	}
	else /* one single and one double solution */
	{
	    double u = cbrt(-q);
	    s[0] = 2 * u;
	    s[1] = - u;
	    num = 2;
	}
    }
    else if (D < 0) /* Casus irreducibilis: three real solutions */
    {
	double phi = double(1.0/3 * acos(float(-q / sqrt(-cb_p))));
	double t = double(2 * sqrt(-p));

	s[0] =   t * cos(float(phi));
	s[1] = double(- t * cos(float(phi + M_PI / 3)));
	s[2] = - t * cos(float(phi - M_PI / 3));
	num = 3;
    }
    else /* one real solution */
    {
	double sqrt_D = sqrt(D);
	double u = cbrt(sqrt_D - q);
	double v = - cbrt(sqrt_D + q);

	s[0] = u + v;
	num = 1;
    }

    /* resubstitute */

    sub = 1.0/3 * A;

    for (i = 0; i < num; ++i)
	s[i] -= sub;

    return num;
}


int SolveQuartic(double c[5], double s[4])
{
    double  coeffs[4];
    double  z, u, v, sub;
    double  A, B, C, D;
    double  sq_A, p, q, r;
    int     i, num;

    /* normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 */

    A = c[3] / c[4];
    B = c[2] / c[4];
    C = c[1] / c[4];
    D = c[0] / c[4];

    /*  substitute x = y - A/4 to eliminate cubic term:
	x^4 + px^2 + qx + r = 0 */

    sq_A = A * A;
    p = - 3.0/8 * sq_A + B;
    q = 1.0/8 * sq_A * A - 1.0/2 * A * B + C;
    r = - 3.0/256*sq_A*sq_A + 1.0/16*sq_A*B - 1.0/4*A*C + D;

    if (IsZero(r))
    {
		/* no absolute term: y(y^3 + py + q) = 0 */

		coeffs[0] = q;
		coeffs[1] = p;
		coeffs[2] = 0;
		coeffs[3] = 1;
	
		double s_cub[3];
		num = SolveCubic(coeffs, s_cub);
		s[0] = s_cub[0];
		s[1] = s_cub[1];
		s[2] = s_cub[2];

		s[num++] = 0;
    }
    else
    {
		/* solve the resolvent cubic ... */

		coeffs[0] = 1.0/2 * r * p - 1.0/8 * q * q;
		coeffs[1] = - r;
		coeffs[2] = - 1.0/2 * p;
		coeffs[3] = 1;

		double s_cub[3];
		num = SolveCubic(coeffs, s_cub);
		s[0] = s_cub[0];
		s[1] = s_cub[1];
		s[2] = s_cub[2];

		/* ... and take the one real solution ... */

		z = s[0];

		/* ... to build two quadric equations */

		u = z * z - r;
		v = 2 * z - p;

		if (IsZero(u))
			u = 0;
		else if (u > 0)
			u = sqrt(u);
		else
			return 0;

		if (IsZero(v))
			v = 0;
		else if (v > 0)
			v = sqrt(v);
		else
			return 0;

		coeffs[0] = z - u;
		coeffs[1] = q < 0 ? -v : v;
		coeffs[2] = 1;

		double s_quad[2];
		double coeff_quad[3];
		coeff_quad[0] = coeffs[0];
		coeff_quad[1] = coeffs[1];
		coeff_quad[2] = coeffs[2];
		num = SolveQuadric(coeff_quad, s_quad);
		int i;
		for(i = 0; i < num; i++)
		{
			s[i] = s_quad[i];
		}

		double s_quad_n[2];
		coeffs[0]= z + u;
		coeffs[1] = q < 0 ? v : -v;
		coeffs[2] = 1;

		coeff_quad[0] = coeffs[0];
		coeff_quad[1] = coeffs[1];
		coeff_quad[2] = coeffs[2];

		num += SolveQuadric(coeff_quad, s_quad_n);
		for(int j = i; j < num; j++)
		{
			s[j] = s_quad[j - i];
		}
    }

    /* resubstitute */

    sub = 1.0/4 * A;

    for (i = 0; i < num; ++i)
	s[i] -= sub;

    return num;
}

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
		vec3 normal = normalize(near - center);

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
CylinderIntersection calcCylinderIntersection(vec3 rayOrigin, vec3 rayDirection, float radius, float height, mat4 transformation)
{
	rayOrigin = vec3(inverse(transformation) * vec4(rayOrigin, 1.0));
	rayDirection = vec3(inverse(transformation) * vec4(rayDirection, 1.0));
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
			return CylinderIntersection(true, vec3(transformation * vec4(rayOrigin + t3 * rayDirection, 1.0f)), vec3( transformation * vec4(near, 1.0f)), vec3(0.0, 1.0, 0.0), vec3(0.0f, 1.0f, 0.0f), t3);
		}
		if (near.y < 0 && 0 < far.y)
		{
		//intersection with lower boundary
			float t3 = (0 - rayOrigin.y) / rayDirection.y;
			return CylinderIntersection(true, vec3(transformation * vec4(rayOrigin + t3 * rayDirection, 1.0f)), vec3( transformation * vec4(near, 1.0f)), vec3(0.0, -1.0, 0.0), vec3(0.0f, -1.0f, 0.0f), t3);
		}
		//calc upper and lower boundary
		if (near.y < 0 || near.y > height)
		{
			return CylinderIntersection(false, near, far, normal_near, normal_far, 0);
		}
		return CylinderIntersection(true, vec3( transformation * vec4(near, 1.0f)), vec3( transformation * vec4(far, 1.0f)), normal_near, normal_far, min(t1,t2));
	}
	else
	{
		return CylinderIntersection(false, vec3(0.0f), vec3(0.0f), vec3(0.0f), vec3(0.0f), 0);
	}
	//check for boundary planes
}

//torus
TorusIntersection calcTorusIntersection(vec3 rayOrigin, vec3 rayDirection, float majR, float minR, mat4 transformation)
{	
	vec3 center = vec3(transformation[3]);
	rayOrigin = vec3(inverse(transformation) * vec4(rayOrigin, 1.0));
	rayDirection = vec3(inverse(transformation) * vec4(rayDirection, 1.0));
	float ddot = dot(rayDirection, rayDirection);
	float odot = dot(rayOrigin, rayOrigin);
	float oddot = dot(rayDirection, rayOrigin);
	double coefficients[5];
	coefficients[4] = double(pow(ddot, 2));
	coefficients[3] = double(4 * ddot * dot(rayOrigin, rayDirection));
	coefficients[2] = double(2 * ddot * (odot - ( (minR * minR) + (majR * majR) )) + 4 * pow(oddot, 2) + 4 * majR * majR * rayDirection.y * rayDirection.y);
	coefficients[1] = double(4 * (odot - ( (minR * minR) + (majR * majR) )) * oddot + 8 * majR * majR * rayOrigin.y * rayDirection.y);
	coefficients[0] = double((odot - ( (minR * minR) + (majR * majR) )) * (odot - ( (minR * minR) + (majR * majR) )) - 4 * majR * majR * (minR * minR - rayOrigin.y * rayOrigin.y));
	//with these parameters we need to solve the equation F(r(t)) = c4 * t^4 + c3 * t^3 + c2 * t^2 + c1 * t + c0 for t and get between 0 and 4 solutions for t
	double solutions[4];
	int num = SolveQuartic(coefficients, solutions);
	if(num > 0) //there is at least one intersection point
	{
		float tMin = maxFloat;
		for(int i = 0; i < num; i++)
		{
			if (solutions[i] < tMin)
			{
				tMin = float(solutions[i]);
			}
		}
		vec3 pos = rayOrigin + tMin * rayDirection;
		vec3 pos_shadow_dir = normalize(vec3(pos.x, 0.0f, pos.z)); //direction from center to the shadow of the point on the xz-plane
		vec3 q =  majR * pos_shadow_dir; //closest point to the point on the center of outer radius
		vec3 normal = normalize(pos - q);
		return TorusIntersection(true, vec3( transformation * vec4(pos, 1.0f)), normal, tMin);
	}
	else
	{
		return TorusIntersection(false, vec3(0.0f), vec3(0.0f), 0);
	}
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
	//vec3 halfway = (light + viewer) / length(light + viewer);
	//float specular = max (dot( halfway, normal), 0.0);
	vec3 reflection = normalize( 2 * dot(light, normalize(normal)) * normalize(normal) - light );
	float specular = max (dot( reflection, viewer), 0.0);

	//calculate illumination
	return vec4( (ka * vec3(0.5, 0.5, 0.5) + kd * dot(light, normalize(normal)) * vec3(0.5f) + ks * pow( specular, shininess) * vec3(1.0f) ), 1.0f);
}
Intersection closestIntersection(Scene scene, vec3 rayOrigin, vec3 rayDirection)
{
	//calculate intersections
	BoxIntersection b = calcBoxIntersection(rayOrigin, rayDirection, scene.box.bmin, scene.box.bmax, scene.box.transformation);
	SphereIntersection s = calcSphereIntersection(scene.sphere.radius, rayOrigin, scene.sphere.center, rayDirection);
	CylinderIntersection c = calcCylinderIntersection(rayOrigin, rayDirection, scene.cylinder.radius, scene.cylinder.height, scene.cylinder.transformation);
	PlaneIntersection p = calcPlaneIntersection(rayOrigin, rayDirection, scene.plane.normal, scene.plane.point);
	//TorusIntersection t = calcTorusIntersection(rayOrigin, rayDirection, scene.torus.major_radius, scene.torus.minor_radius, scene.torus.transformation);

	//find closest intersection -> minimal t_near value, if there is one
	float tMin = maxFloat;
	int shininess;
	vec3 normal, ka, kd, ks;
	vec4 color;
	if(! (s.hit || c.hit || p.hit || b.hit /*|| t.hit*/))
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
		shininess = scene.sphere.shininess;
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
	if(p.hit && p.t_near < tMin)
	{
		tMin = p.t_near;
		normal = scene.plane.normal;
		color = vec4(0.7f, 0.7f, 0.7f, 1.0f);
		ka = scene.plane.ka;
		kd = scene.plane.kd;
		ks = scene.plane.ks;
		shininess = scene.plane.shininess;
	}
	/*if(t.hit && t.t_near < tMin)
	{
		tMin = t.t_near;
		normal = t.normal_near;
		color = vec4(0.0f, 0.7f, 0.0f, 1.0f);
		ka = scene.torus.ka;
		kd = scene.torus.kd;
		ks = scene.torus.ks;
		shininess = scene.torus.shininess;
	}*/

	return Intersection(true, tMin, normal, color, ka, kd, ks, shininess);
}

void main()
{
	const float _zero=1e-6;
	vec4 near = inverseModelViewProjectionMatrix*vec4(fragPosition,-1.0,1.0); //world space ray position (from screen space back to world space)
	near /= near.w;	//perspective divide -> diverge rays

	vec4 far = inverseModelViewProjectionMatrix*vec4(fragPosition,1.0,1.0);
	far /= far.w;

	// this is the setup for our viewing ray
	vec3 rayOrigin = near.xyz;
	vec3 rayDirection = normalize((far-near).xyz);

	float idx = (sin(0.5 * currentTime) + 1) / 2.0f;
	vec3 spherePos = cubicBezierVector(vec3(1.0f), vec3(1.0f, -1.0f, 0.0f), vec3(1.0f, -1.0f, -1.0f), vec3 (1.0f), idx);
	vec3 boxPos = cubicBezierVector(vec3(-1.0f), vec3(-0.5f, 0.5f, 0.0f), vec3(-2.0f, 1.0f, -1.5f), vec3 (-1.0f), idx);

	mat4 boxTransformation = mat4(1.0f);
	boxTransformation[0][0] = cos(radians(45));
	boxTransformation[1][1] = cos(radians(45));
	boxTransformation[1][0] = - sin(radians(45));
	boxTransformation[0][1] = sin(radians(45));

	mat4 cylinderTransformation = mat4(1.0f);
	cylinderTransformation[0][0] = cos(radians(30));
	cylinderTransformation[1][1] = cos(radians(30));
	cylinderTransformation[1][0] = - sin(radians(30));
	cylinderTransformation[0][1] = sin(radians(30));
	cylinderTransformation[3] = vec4(-0.1f, 0.05f, 0.15f, 1.0f);
	vec3 ka_mtl = vec3(0.0f);
	vec3 kd_mtl = vec3(0.7f);
	vec3 ks_mtl = vec3(0.2f);
	int shininess = 50;

	Box box = Box(boxPos, boxPos + 0.5, boxTransformation, ka_mtl, kd_mtl, ks_mtl, shininess);
	Sphere sphere = Sphere(spherePos, 0.5, ka_mtl, kd_mtl, ks_mtl, shininess);
	Plane plane = Plane(vec3(0.0f, 1.0f, 0.0f), vec3(-1.0f, -1.0f, -1.0f), ka_mtl, kd_mtl, ks_mtl, shininess);
	Cylinder cylinder = Cylinder(0.4, 0.7, cylinderTransformation, ka_mtl, kd_mtl, ks_mtl, shininess);
	//Torus torus = Torus(0.6f, 0.2f, mat4(1.0f), ka_mtl, kd_mtl, ks_mtl, shininess);

	Scene scene = Scene(box, sphere, plane, cylinder);
	
	if (enableRT)
	{
		Intersection sec = closestIntersection(scene, rayOrigin, rayDirection);
		if(sec.hit)
		{
			vec3 n = sec.normal;
			vec3 p = rayOrigin + sec.tMin * rayDirection;
			fragColor = sec.color * phongIllumination(p, rayOrigin, n, sec.ka, sec.kd, sec.ks, sec.shininess);
			gl_FragDepth = calcDepth(p);
			return;	
		}
	}
	
	//raytracing setup -> commented out for performance and testing of scene
	/*
	const int max_bounces = 4;
	const int max_rays = (1 << max_bounces) - 1;
	ray rays[max_rays];
	int num_rays;
	vec3 position, light, v0;
	float t, n1;

	//start ray in array
	rays[0].p = rayOrigin;
	rays[0].d = normalize(rayDirection);
	rays[0].n = vec3(0.0f);
	rays[0].rfl = 0.0f;
	rays[0].rfr = 0.0f;
	rays[0].n0 = n0;
	rays[0].n1 = 1.54f;
	rays[0].l = 0.0f;
	rays[0].level = 0;
	rays[0].i0 = -1;
	rays[0].i1 = -1;
	rays[0].sec = Intersection(false, 0, vec3(0), vec4(0), vec3(0), vec3(0), vec3(0), 0);
	num_rays = 1;

	for( int i = 0; i < num_rays; i++)
	{
		//get the closest intersection and its properties
		n1 = 1.54f;
		Intersection sec = closestIntersection(scene, rays[i].p, rays[i].d);
		rays[i].sec = sec;
		rays[i].color = sec.color; //color of intersected volume
		rays[i].rfl = kr; //reflection value
		rays[i].rfr = kt; //refraction value
		rays[i].n = normalize(sec.normal); //normal at intersection
		rays[i].l = sec.tMin; //length of ray
		position = rays[i].p + rays[i].l * rays[i].d;

		//reflect and refract ray
		if (sec.hit && rays[i].level < max_bounces - 1)
		{
			//reflecion
			t = dot(rays[i].d, rays[i].n);
			if (rays[i].rfl > _zero && t < _zero)
			{
				light = normalize( worldLightPosition - position );
				rays[i].i0 = num_rays;
				rays[num_rays] = rays[i];
				rays[num_rays].level++;
				rays[num_rays].i0 = -1;
				rays[num_rays].i1 = -1;
				rays[num_rays].p = position;
				rays[num_rays].d = normalize( 2 * dot(light, normalize(rays[i].n)) * normalize(rays[i].n) - light );
				rays[num_rays].n0 = rays[i].n0;
				rays[num_rays].n1 = rays[0].n0;
				num_rays++;
			}
			//refraction
			if (rays[i].rfr > _zero)
			{
				rays[i].i1 = num_rays;
				rays[num_rays] = rays[i];
				rays[num_rays].level++;
				rays[num_rays].i0 = -1;
				rays[num_rays].i1 = -1;
				rays[num_rays].p = position;
				if (t > 0) //exit the object
				{
					rays[num_rays].n0 = rays[i].n0;
					rays[num_rays].n1 = n0;
					v0 = -rays[i].n;
					t = -t;
				}
				else //enter the object
				{
					rays[num_rays].n0 = n1;
					rays[num_rays].n1 = rays[i].n0;
					rays[i].n1 = n1;
					v0 = rays[i].n;
				}
				n1 = rays[i].n0 / rays[i].n1;
				t = 1.0f - (n1 * n1 * (1.0f - t * t));
				if (t > 0)
				{
					rays[num_rays].d = (rays[i].d * n1) - (v0 * ((n1 * t) + sqrt(t)));
				}
				num_rays++;
			}
		}
		else //ignore ray if nothing hit
		{
			rays[i] = rays[num_rays-1];
			num_rays--;
			i--;
		}
	}
	//backtrack of the intersections
	for (int i = num_rays - 1; i >= 0; i--)
	{
		// directional + ambient light
		t = phongIllumination(rays[i].p + rays[i].l * rays[i].d, rays[i].p, rays[i].n, rays[i].sec.ka, rays[i].sec.kd, rays[i].sec.ks, rays[i].sec.shininess);
		t *= 1.0-rays[i].rfl-rays[i].rfr;
		rays[i].color.rgb *= t;
		// reflect
		int _i=rays[i].i0;
		if (_i >= 0) 
		{
			rays[i].color.rgb += rays[_i].color.rgb * rays[i].rfl;
		}
		// refract
		_i = rays[_i].i1;
		if (_i >= 0) 
		{
			rays[i].color.rgb += rays[_i].color.rgb * rays[i].rfr;
		}
	}
	fragColor = rays[0].color;
	gl_FragDepth = calcDepth(rays[0].p + rays[0].l * rays[0].d);
	return;
	*/
	// using calcDepth, you can convert a ray position to an OpenGL z-value, so that intersections/occlusions with the
	// model geometry are handled correctly, e.g.: gl_FragDepth = calcDepth(nearestHit);
	// in case there is no intersection, you should get gl_FragDepth to 1.0, i.e., the output of the shader will be ignored
	fragColor = vec4(1.0);
	gl_FragDepth = 1.0; //-> visibility and intersection
}
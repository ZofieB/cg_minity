#version 400
#extension GL_ARB_shading_language_include : require
#include "/model-globals.glsl"

uniform vec3 worldCameraPosition;
uniform vec3 worldLightPosition;
uniform sampler2D diffuseTexture;
uniform bool wireframeEnabled;
uniform vec4 wireframeLineColor;

//new uniforms
uniform vec3 ks = vec3(1.0f, 1.0f, 1.0f);
uniform vec3 ka = vec3(1.0f, 1.0f, 1.0f);
uniform vec3 kd = vec3(1.0f, 1.0f, 1.0f);
uniform vec3 ia;
uniform vec3 is;
uniform vec3 id;
uniform int shininess;
uniform mat3 normalMatrix;


in fragmentData
{
	vec3 position;
	vec3 normal;
	vec2 texCoord;
	noperspective vec3 edgeDistance;
} fragment;

out vec4 fragColor;

void main()
{
	//transform fragment normal 
	vec3 normal = normalize(fragment.normal);

	//direction vector from surface point to light source
	vec3 light = normalize( worldLightPosition - fragment.position );

	//direction vector pointing towards viewer
	vec3 viewer = normalize( worldCameraPosition - fragment.position );

	//calculate reflecition vector, reflection only if in viewer direction (positive), else it is 0
	vec3 reflection = normalize( 2 * dot(light, normal) * normal - light );
	float specular = max (dot( reflection, viewer), 0.0);

	//calculate result of phong illumination
	vec4 result = vec4( (ka * ia + kd * dot(light, normal) * id + ks * pow( specular, shininess) * is ), 1.0f);

	if (wireframeEnabled)
	{
		float smallestDistance = min(min(fragment.edgeDistance[0],fragment.edgeDistance[1]),fragment.edgeDistance[2]);
		float edgeIntensity = exp2(-1.0*smallestDistance*smallestDistance);
		result.rgb = mix(result.rgb,wireframeLineColor.rgb,edgeIntensity*wireframeLineColor.a);
	}

	fragColor = result;
}
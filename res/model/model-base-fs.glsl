#version 400
#extension GL_ARB_shading_language_include : require
#include "/model-globals.glsl"

uniform vec3 worldCameraPosition;
uniform vec3 worldLightPosition;
uniform bool wireframeEnabled;
uniform vec4 wireframeLineColor;

//new uniforms -- assignment 1
uniform vec3 ks = vec3(1.0f, 1.0f, 1.0f);
uniform vec3 ka = vec3(1.0f, 1.0f, 1.0f);
uniform vec3 kd = vec3(1.0f, 1.0f, 1.0f);
uniform vec3 ia;
uniform vec3 is;
uniform vec3 id;
uniform int gui_shininess;
uniform float mat_shininess = - 1.0;
uniform mat3 normalMatrix;

//new variables -- assignment 2
uniform sampler2D diffuseTexture;
uniform bool diff_tex_loaded = false;
uniform sampler2D ambientTexture;
uniform bool amb_tex_loaded = false;
uniform sampler2D specularTexture;
uniform bool spec_tex_loaded = false;

uniform sampler2D objectSpaceNormalTexture;
uniform sampler2D tangentSpaceNormalTexture;

uniform bool diffuseTextureEnabled;
uniform bool ambientTextureEnabled;
uniform bool specularTextureEnabled;
uniform bool obj_norm_map;
uniform bool tan_norm_map;

uniform mat4 modelViewProjectionMatrix;

in fragmentData
{
	vec3 position;
	vec3 normal;
	vec2 texCoord;
	vec3 tangent;
	vec3 bitangent;
	noperspective vec3 edgeDistance;
} fragment;

out vec4 fragColor;

void main()
{
	//compute fragment normal -> normal mapping
	vec3 normal;
	if (obj_norm_map)
	{
		normal = texture(objectSpaceNormalTexture, fragment.texCoord).rgb;
		normal = normalize(normal * 2.0 - 1.0);
	}
	else if (tan_norm_map)
	{
		normal = texture(tangentSpaceNormalTexture, fragment.texCoord).rgb;
		normal = normal * 2.0 - 1.0;

		vec3 T = normalize( normalMatrix * vec3(modelViewProjectionMatrix * vec4(fragment.tangent, 0.0)));
		//vec3 B = normalize( normalMatrix * vec3(modelViewProjectionMatrix * vec4(fragment.bitangent, 0.0)));
		vec3 N = normalize(normal);

		//using Gram-Schmidt to re-orthogonalize
		T = normalize(T - dot(T, N) * N);
		vec3 B = cross(N, T);

		mat3 TBN = mat3(T, B, N);
		normal = normalize(TBN * normal);
	}
	else
	{
		normal = normalize(fragment.normal);
	}

	//if material shininess is loaded
	float shininess;
	if (mat_shininess != - 1.0)
	{
		shininess = mat_shininess;
	}
	else
	{
		shininess = gui_shininess;
	}
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
	// TODO: dont load textures, if they are not specified! (even if checkbox active)
	if (diffuseTextureEnabled && diff_tex_loaded)
	{
		 vec4 diffuse_tex = texture(diffuseTexture, fragment.texCoord);
		 result = result * diffuse_tex;
	}
	if (ambientTextureEnabled && amb_tex_loaded)
	{
		vec4 ambient_tex = texture(ambientTexture, fragment.texCoord);
		result = result * ambient_tex;
	}
	if (specularTextureEnabled && spec_tex_loaded)
	{
		vec4 specular_tex = texture(specularTexture, fragment.texCoord);
		result = result * specular_tex;
	}

	fragColor = result;
}
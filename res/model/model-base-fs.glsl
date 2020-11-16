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
uniform bool diffuseTextureEnabled;

uniform sampler2D ambientTexture;
uniform bool amb_tex_loaded = false;
uniform bool ambientTextureEnabled;

uniform sampler2D specularTexture;
uniform bool spec_tex_loaded = false;
uniform bool specularTextureEnabled;

uniform sampler2D objectSpaceNormalTexture;
uniform bool obj_norm_map;
uniform bool obj_norm_map_loaded = false;
uniform sampler2D tangentSpaceNormalTexture;
uniform bool tan_norm_map;
uniform bool tan_norm_map_loaded = false;

uniform bool bump_map;
uniform bool sinus;
uniform bool sinus_sqr;
uniform float a;
uniform float k;

in fragmentData
{
	vec3 position;
	vec3 normal;
	vec2 texCoord;
	mat3 TBNMatrix;
	vec3 tangent;
	vec3 bitangent;
	noperspective vec3 edgeDistance;
} fragment;

out vec4 fragColor;

void main()
{
	//compute fragment normal -> normal mapping
	vec3 normal;
	vec3 position = fragment.position;
	//object space normal mapping
	if (obj_norm_map_loaded && obj_norm_map)
	{
		normal = texture(objectSpaceNormalTexture, fragment.texCoord).rgb;
		normal = normalize(normal * 2.0 - 1.0);
	}
	//tangent space normal mapping
	else if (tan_norm_map_loaded && tan_norm_map)
	{
		normal = texture(tangentSpaceNormalTexture, fragment.texCoord).rgb;
		normal = normalize(normal * 2.0 - 1.0);
		if (bump_map)
		{
			float bump_map_value, bu, bv;
			if (sinus_sqr)
			{
				//bump map function
				bump_map_value = a * (sin(k * fragment.texCoord.x) * sin(k * fragment.texCoord.x) * sin(k * fragment.texCoord.y) * sin(k * fragment.texCoord.y));
				//derivative bump map in u direction
				bu = dFdx(bump_map_value);
				//derivative bump map in v direction
				bv = dFdy(bump_map_value);
			}
			else if (sinus)
			{
				//other bump map function
				bump_map_value = a * sin(k * fragment.texCoord.x) * sin(k * fragment.texCoord.y);
				bu = a * sin(k * fragment.texCoord.y) * cos(k * fragment.texCoord.x) * k;
				bv = a * sin(k * fragment.texCoord.x) * cos(k * fragment.texCoord.y) * k;
			}
			//transform normal according to bump map
			normal = normalize(normal + bv * cross(fragment.tangent, normal) + bu * cross(normal, fragment.bitangent));
		}
		//transform normal from tangent space to world space
		normal = normalize(fragment.TBNMatrix * normal);
	}
	else
	{
		normal = normalize(fragment.normal);
	}

	//pick correct shininess based on given material shininess
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
	vec3 light = normalize( worldLightPosition - position );

	//direction vector pointing towards viewer
	vec3 viewer = normalize( worldCameraPosition - position );

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
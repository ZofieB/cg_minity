#version 400
#extension GL_ARB_shading_language_include : require
#include "/model-globals.glsl"

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

uniform vec2 viewportSize;
uniform mat3 normalMatrix;

in vertexData
{
	vec3 position;
	vec3 normal;
	vec2 texCoord;
} vertices[];

out fragmentData
{
	vec3 position;
	vec3 normal;
	vec2 texCoord;
	mat3 TBNMatrix;
	vec3 tangent;
	vec3 bitangent;
	noperspective vec3 edgeDistance;
} fragment;

void main(void)
{
	vec2 p[3];
	vec2 v[3];

	for (int i=0;i<3;i++)
		p[i] = 0.5 * viewportSize *  gl_in[i].gl_Position.xy/gl_in[i].gl_Position.w;

	v[0] = p[2]-p[1];
	v[1] = p[2]-p[0];
	v[2] = p[1]-p[0];

	float area = abs(v[1].x*v[2].y - v[1].y * v[2].x);

	for (int i=0;i<3;i++)
	{
		gl_Position = gl_in[i].gl_Position;
		fragment.position = vertices[i].position;
		fragment.normal = vertices[i].normal;
		fragment.texCoord = vertices[i].texCoord;
		
		vec3 ed = vec3(0.0);
		ed[i] = area / length(v[i]);
		fragment.edgeDistance = ed;

		
		//calculate TBNMatrix
		//triangle edges and uv distances 
		vec3 edge1 = vertices[(i + 1) % 3].position - vertices[i].position;
		vec3 edge2 = vertices[(i + 2) % 3].position - vertices[i].position;
		vec2 deltaUV1 = vertices[(i + 1) % 3].texCoord - vertices[i].texCoord;
		vec2 deltaUV2 = vertices[(i + 2) % 3].texCoord - vertices[i].texCoord;

		//calculate tangent and bitangent out of edges and texture coordinates

		vec3 tangent, bitangent;

		float invDet = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

		tangent.x = invDet * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
		tangent.y = invDet * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
		tangent.z = invDet * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);

		bitangent.x = invDet * (- deltaUV2.x * edge1.x + deltaUV1.x * edge2.x);
		bitangent.y = invDet * (- deltaUV2.x * edge1.y + deltaUV1.x * edge2.y);
		bitangent.z = invDet * (- deltaUV2.x * edge1.z + deltaUV1.x * edge2.z);

		tangent = normalize(normalMatrix * tangent);
		bitangent = normalize(normalMatrix * bitangent);
		vec3 normal = normalize(normalMatrix * vertices[i].normal);

		fragment.TBNMatrix = mat3(tangent, bitangent, normal);
		fragment.tangent = tangent;
		fragment.bitangent = bitangent;

		EmitVertex();
	}

	EndPrimitive();
}
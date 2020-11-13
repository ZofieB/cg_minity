#include "ModelRenderer.h"
#include <globjects/base/File.h>
#include <globjects/State.h>
#include <iostream>
#include <filesystem>
#include <imgui.h>
#include "Viewer.h"
#include "Scene.h"
#include "Model.h"
#include <sstream>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

using namespace minity;
using namespace gl;
using namespace glm;
using namespace globjects;

ModelRenderer::ModelRenderer(Viewer* viewer) : Renderer(viewer)
{
	m_lightVertices->setStorage(std::array<vec3, 1>({ vec3(0.0f) }), GL_NONE_BIT);
	auto lightVertexBinding = m_lightArray->binding(0);
	lightVertexBinding->setBuffer(m_lightVertices.get(), 0, sizeof(vec3));
	lightVertexBinding->setFormat(3, GL_FLOAT);
	m_lightArray->enable(0);
	m_lightArray->unbind();

	createShaderProgram("model-base", {
		{ GL_VERTEX_SHADER,"./res/model/model-base-vs.glsl" },
		{ GL_GEOMETRY_SHADER,"./res/model/model-base-gs.glsl" },
		{ GL_FRAGMENT_SHADER,"./res/model/model-base-fs.glsl" },
		}, 
		{ "./res/model/model-globals.glsl" });

	createShaderProgram("model-light", {
		{ GL_VERTEX_SHADER,"./res/model/model-light-vs.glsl" },
		{ GL_FRAGMENT_SHADER,"./res/model/model-light-fs.glsl" },
		}, { "./res/model/model-globals.glsl" });
}

void ModelRenderer::display()
{
	// Save OpenGL state
	auto currentState = State::currentState();

	// retrieve/compute all necessary matrices and related properties
	const mat4 viewMatrix = viewer()->viewTransform();
	const mat4 inverseViewMatrix = inverse(viewMatrix);
	const mat4 modelViewMatrix = viewer()->modelViewTransform();
	const mat4 inverseModelViewMatrix = inverse(modelViewMatrix);
	const mat4 modelLightMatrix = viewer()->modelLightTransform();
	const mat4 inverseModelLightMatrix = inverse(modelLightMatrix);
	const mat4 modelViewProjectionMatrix = viewer()->modelViewProjectionTransform();
	const mat4 inverseModelViewProjectionMatrix = inverse(modelViewProjectionMatrix);
	const mat4 projectionMatrix = viewer()->projectionTransform();
	const mat4 inverseProjectionMatrix = inverse(projectionMatrix);
	const mat3 normalMatrix = mat3(transpose(inverseModelViewMatrix));
	const mat3 inverseNormalMatrix = inverse(normalMatrix);
	const vec2 viewportSize = viewer()->viewportSize();

	auto shaderProgramModelBase = shaderProgram("model-base");

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	viewer()->scene()->model()->vertexArray().bind();

	const std::vector<Group> & groups = viewer()->scene()->model()->groups();
	const std::vector<Material> & materials = viewer()->scene()->model()->materials();

	static std::vector<bool> groupEnabled(groups.size(), true);
	static bool wireframeEnabled = false;
	static bool lightSourceEnabled = true;
	static vec4 wireframeLineColor = vec4(1.0f);
	//static bool modelEnabled = true;

	if (ImGui::BeginMenu("Model"))
	{
		ImGui::Checkbox("Wireframe Enabled", &wireframeEnabled);
		ImGui::Checkbox("Light Source Enabled", &lightSourceEnabled);
		ImGui::Checkbox("Model Enabled", &(viewer()->modelEnabled));

		if (wireframeEnabled)
		{
			if (ImGui::CollapsingHeader("Wireframe"))
			{
				ImGui::ColorEdit4("Line Color", (float*)&wireframeLineColor, ImGuiColorEditFlags_AlphaBar);
			}
		}

		if (ImGui::CollapsingHeader("Groups"))
		{
			for (uint i = 0; i < groups.size(); i++)
			{
				bool checked = groupEnabled.at(i);
				ImGui::Checkbox(groups.at(i).name.c_str(), &checked);
				groupEnabled[i] = checked;
			}

		}

		ImGui::EndMenu();
	}

	if (ImGui::Begin("Parameter"))
	{
		ImGui::Text("Assignment 1");
		ImGui::SliderInt("shininess", &m_shininess, 5, 100);
		ImGui::ColorEdit3("ia", (float*)&m_ambient);
		ImGui::ColorEdit3("is", (float*)&m_specular);
		ImGui::ColorEdit3("id", (float*)&m_diffuse);

		ImGui::Text("Assignment 2");
		ImGui::Checkbox("Diffuse texturing", &diff_tex);
		ImGui::Checkbox("Ambient texturing", &amb_tex);
		ImGui::Checkbox("Specular texturing", &spec_tex);

		if (ImGui::BeginMenu("Normal Maps"))
		{
			if (ImGui::MenuItem("None"))
			{
				obj_norm_map = false;
				tan_norm_map = false;
			}
			if (ImGui::MenuItem("Object Space Normal Mapping"))
			{
				obj_norm_map = true;
				tan_norm_map = false;
			}
			if (ImGui::MenuItem("Tangent Space Normal Mapping"))
			{
				tan_norm_map = true;
				obj_norm_map = false;
			}
			ImGui::EndMenu();
		}
		if (tan_norm_map)
		{
			ImGui::Checkbox("Bump mapping", &bump_map);
			if (bump_map)
			{
				if (ImGui::BeginMenu("Bump map function"))
				{
					if (ImGui::MenuItem("Sinus function"))
					{
						sinus = true;
						sinus_sqr = false;
					}
					if (ImGui::MenuItem("Sinus squared function"))
					{
						sinus_sqr = true;
						sinus = false;
					}
					ImGui::EndMenu();
				}
				if (sinus || sinus_sqr)
				{
					ImGui::SliderFloat("A", &amplitude, 0.1f, 8.0f);
					ImGui::SliderFloat("k", &frequency, 5.0f, 500.0f);
				}
			}
		}

		ImGui::Text("Assignment 3");
		ImGui::SliderFloat("Explosion factor", &viewer()->m_explosion, 0.0f, 50.0f);

		ImGui::End();
	}

	vec4 worldCameraPosition = inverseModelViewMatrix * vec4(0.0f, 0.0f, 0.0f, 1.0f);
	vec4 worldLightPosition = inverseModelLightMatrix * vec4(0.0f, 0.0f, 0.0f, 1.0f);

	shaderProgramModelBase->setUniform("modelViewProjectionMatrix", modelViewProjectionMatrix);
	shaderProgramModelBase->setUniform("modelViewMatrix", modelViewMatrix);
	shaderProgramModelBase->setUniform("viewportSize", viewportSize);
	shaderProgramModelBase->setUniform("worldCameraPosition", vec3(worldCameraPosition));
	shaderProgramModelBase->setUniform("worldLightPosition", vec3(worldLightPosition));
	shaderProgramModelBase->setUniform("wireframeEnabled", wireframeEnabled);
	shaderProgramModelBase->setUniform("wireframeLineColor", wireframeLineColor);

	//illumination model
	shaderProgramModelBase->setUniform("normalMatrix", normalMatrix);
	shaderProgramModelBase->setUniform("gui_shininess", m_shininess);
	shaderProgramModelBase->setUniform("ia", m_ambient);
	shaderProgramModelBase->setUniform("is", m_specular);
	shaderProgramModelBase->setUniform("id", m_diffuse);

	//bump mapping
	shaderProgramModelBase->setUniform("k", frequency);
	shaderProgramModelBase->setUniform("a", amplitude);

	//exploded view
	shaderProgramModelBase->setUniform("explosion_factor", viewer()->m_explosion);

	//boolean flags for textures, bump maps and GUI control
	shaderProgramModelBase->setUniform("diffuseTextureEnabled", diff_tex);
	shaderProgramModelBase->setUniform("ambientTextureEnabled", amb_tex);
	shaderProgramModelBase->setUniform("specularTextureEnabled", spec_tex);
	shaderProgramModelBase->setUniform("obj_norm_map", obj_norm_map);
	shaderProgramModelBase->setUniform("tan_norm_map", tan_norm_map);

	shaderProgramModelBase->setUniform("bump_map", bump_map);
	shaderProgramModelBase->setUniform("sinus", sinus);
	shaderProgramModelBase->setUniform("sinus_sqr", sinus_sqr);
	
	shaderProgramModelBase->use();

	for (uint i = 0; i < groups.size(); i++)
	{
		//use i index to access offset and pass to shader
		if (groupEnabled.at(i))
		{
			const Material & material = materials.at(groups.at(i).materialIndex);
			const Group& group = groups.at(i);

			shaderProgramModelBase->setUniform("group_center", group.center);

			//pass material properties to shader for phong illumintation
			shaderProgramModelBase->setUniform("kd", material.diffuse);
			shaderProgramModelBase->setUniform("ks", material.specular);
			shaderProgramModelBase->setUniform("ka", material.ambient);
			//load material shininess if present instead of shininess slider
			if (material.shininess)
			{
				shaderProgramModelBase->setUniform("mat_shininess", material.shininess);
			}

			if (material.diffuseTexture)
			{
				shaderProgramModelBase->setUniform("diff_tex_loaded", true);
				shaderProgramModelBase->setUniform("diffuseTexture", 0);
				material.diffuseTexture->bindActive(0);
			}
			if (material.ambientTexture)
			{
				shaderProgramModelBase->setUniform("amb_tex_loaded", true);
				shaderProgramModelBase->setUniform("ambientTexture", 1);
				material.ambientTexture->bindActive(1);
			}
			if (material.specularTexture)
			{
				shaderProgramModelBase->setUniform("spec_tex_loaded", true);
				shaderProgramModelBase->setUniform("specularTexture", 2);
				material.specularTexture->bindActive(2);
			}
			if (material.objectSpaceNormalTexture)
			{
				shaderProgramModelBase->setUniform("obj_norm_map_loaded", true);
				shaderProgramModelBase->setUniform("objectSpaceNormalTexture", 3);
				material.objectSpaceNormalTexture->bindActive(3);
			}
			if (material.tangentSpaceNormalTexture)
			{
				shaderProgramModelBase->setUniform("tan_norm_map_loaded", true);
				shaderProgramModelBase->setUniform("tangentSpaceNormalTexture", 4);
				material.tangentSpaceNormalTexture->bindActive(4);
			}

			if (viewer()->modelEnabled)
			{
				viewer()->scene()->model()->vertexArray().drawElements(GL_TRIANGLES, groups.at(i).count(), GL_UNSIGNED_INT, (void*)(sizeof(GLuint) * groups.at(i).startIndex));
			}

			if (material.tangentSpaceNormalTexture)
			{
				material.tangentSpaceNormalTexture->unbind();
			}
			if (material.objectSpaceNormalTexture)
			{
				material.objectSpaceNormalTexture->unbind();
			}
			if (material.specularTexture)
			{
				material.specularTexture->unbind();
			}
			if (material.ambientTexture)
			{
				material.ambientTexture->unbind();
			}
			if (material.diffuseTexture)
			{
				material.diffuseTexture->unbind();
			}

		}
	}

	shaderProgramModelBase->release();

	viewer()->scene()->model()->vertexArray().unbind();


	if (lightSourceEnabled)
	{
		auto shaderProgramModelLight = shaderProgram("model-light");
		shaderProgramModelLight->setUniform("modelViewProjectionMatrix", modelViewProjectionMatrix * inverseModelLightMatrix);
		shaderProgramModelLight->setUniform("viewportSize", viewportSize);

		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDepthMask(GL_FALSE);

		m_lightArray->bind();

		shaderProgramModelLight->use();
		m_lightArray->drawArrays(GL_POINTS, 0, 1);
		shaderProgramModelLight->release();

		m_lightArray->unbind();

		glDisable(GL_PROGRAM_POINT_SIZE);
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
	}

	// Restore OpenGL state (disabled to to issues with some Intel drivers)
	// currentState->apply();
}
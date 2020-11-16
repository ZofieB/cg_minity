#include "RaytraceRenderer.h"
#include <globjects/base/File.h>
#include <globjects/State.h>
#include <iostream>
#include <filesystem>
#include <imgui.h>
#include "Viewer.h"
#include "Scene.h"
#include "Model.h"
#include <sstream>
#include <limits>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

using namespace minity;
using namespace gl;
using namespace glm;
using namespace globjects;

RaytraceRenderer::RaytraceRenderer(Viewer* viewer) : Renderer(viewer)
{
	m_quadVertices->setStorage(std::array<vec2, 4>({ vec2(-1.0f, 1.0f), vec2(-1.0f,-1.0f), vec2(1.0f,1.0f), vec2(1.0f,-1.0f) }), gl::GL_NONE_BIT);
	auto vertexBindingQuad = m_quadArray->binding(0);
	vertexBindingQuad->setBuffer(m_quadVertices.get(), 0, sizeof(vec2));
	vertexBindingQuad->setFormat(2, GL_FLOAT);
	m_quadArray->enable(0);
	m_quadArray->unbind();

	createShaderProgram("raytrace", {
			{ GL_VERTEX_SHADER,"./res/raytrace/raytrace-vs.glsl" },
			{ GL_FRAGMENT_SHADER,"./res/raytrace/raytrace-fs.glsl" },
		}, 
		{ "./res/raytrace/raytrace-globals.glsl" });
}

void RaytraceRenderer::display()
{
	// Save OpenGL state
	auto currentState = State::currentState();

	// retrieve/compute all necessary matrices and related properties
	const mat4 modelViewProjectionMatrix = viewer()->modelViewProjectionTransform();
	const mat4 inverseModelViewProjectionMatrix = inverse(modelViewProjectionMatrix);
	const mat4 modelViewMatrix = viewer()->modelViewTransform();
	const mat4 inverseModelViewMatrix = inverse(modelViewMatrix);
	const mat4 modelLightMatrix = viewer()->modelLightTransform();
	const mat4 inverseModelLightMatrix = inverse(modelLightMatrix);

	vec4 worldLightPosition = inverseModelLightMatrix * vec4(0.0f, 0.0f, 0.0f, 1.0f);
	vec4 worldCameraPosition = inverseModelViewMatrix * vec4(0.0f, 0.0f, 0.0f, 1.0f);

	float currentTime = (float)glfwGetTime();

	auto shaderProgramRaytrace = shaderProgram("raytrace");

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	float mx_kr = 0.3;
	float max_kl = 0.4;
	float max_kt = 0.3;

	if (ImGui::Begin("Parameter Raytracer"))
	{
		ImGui::Checkbox("Enable RT-scene", &raytracing_enabled);
		ImGui::SliderFloat("Reflection constant", &kr, 0.0f, 1.0f - kl - kt);
		ImGui::SliderFloat("Illumination constant", &kl, 0.0f, 1.0f - kr - kt);
		ImGui::SliderFloat("Refraction constant", &kt, 0.0f, 1.0f - kr - kl);
		ImGui::End();
	}

	shaderProgramRaytrace->setUniform("modelViewProjectionMatrix", modelViewProjectionMatrix);
	shaderProgramRaytrace->setUniform("inverseModelViewProjectionMatrix", inverseModelViewProjectionMatrix);
	shaderProgramRaytrace->setUniform("worldCameraPosition", vec3(worldCameraPosition));
	shaderProgramRaytrace->setUniform("worldLightPosition", vec3(worldLightPosition));
	shaderProgramRaytrace->setUniform("currentTime", currentTime);
	shaderProgramRaytrace->setUniform("maxFloat", std::numeric_limits<float>::max());
	shaderProgramRaytrace->setUniform("kr", kr);
	shaderProgramRaytrace->setUniform("kl", kl);
	shaderProgramRaytrace->setUniform("kt", kt);
	shaderProgramRaytrace->setUniform("enableRT", raytracing_enabled);

	m_quadArray->bind();
	shaderProgramRaytrace->use();
	// we are rendering a screen filling quad (as a tringle strip), so we can cast rays for every pixel
	m_quadArray->drawArrays(GL_TRIANGLE_STRIP, 0, 4);
	shaderProgramRaytrace->release();
	m_quadArray->unbind();


	// Restore OpenGL state (disabled to to issues with some Intel drivers)
	// currentState->apply();
}
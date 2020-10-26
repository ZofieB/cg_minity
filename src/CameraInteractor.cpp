#include "CameraInteractor.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <time.h>
#include <chrono>
#include <thread>

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glbinding/gl/gl.h>
#include <glbinding/gl/enum.h>
#include <glbinding/gl/functions.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include "Viewer.h"
#include "ModelRenderer.h"

using namespace minity;
using namespace glm;
using namespace gl;

CameraInteractor::CameraInteractor(Viewer * viewer) : Interactor(viewer)
{
	resetProjectionTransform();
	resetViewTransform();

	globjects::debug() << "Camera interactor usage:";
	globjects::debug() << "  Drag left mouse - rotate";
	globjects::debug() << "  Drag middle mouse - pan";
	globjects::debug() << "  Drag right mouse - zoom";
	globjects::debug() << "  Shift + Left mouse - light position";
	globjects::debug() << "  H - toggle headlight";
	globjects::debug() << "  B - benchmark";
	globjects::debug() << "  Home - reset view";
	globjects::debug() << "  Cursor left - rotate negative around current y-axis";
	globjects::debug() << "  Cursor right - rotate positive around current y-axis";
	globjects::debug() << "  Cursor up - rotate negative around current x-axis";
	globjects::debug() << "  Cursor right - rotate positive around current x-axis" << std::endl;
}

void CameraInteractor::framebufferSizeEvent(int width, int height)
{
	float aspect = float(width) / float(height);

	if (m_perspective)
		viewer()->setProjectionTransform(perspective(m_fov, aspect, m_near, m_far));
	else
		viewer()->setProjectionTransform(ortho(-1.0f*aspect, 1.0f*aspect, -1.0f, 1.0f, m_near, m_far));
}

void CameraInteractor::keyEvent(int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS)
	{
		m_light = true;
		m_xPrevious = m_xCurrent;
		m_yPrevious = m_yCurrent;
		cursorPosEvent(m_xCurrent, m_yCurrent);
	}
	else if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE)
	{
		m_light = false;
	}
	else if (key == GLFW_KEY_HOME && action == GLFW_RELEASE)
	{
		resetViewTransform();
	}
	else if (key == GLFW_KEY_LEFT && action == GLFW_RELEASE)
	{
		mat4 viewTransform = viewer()->viewTransform();
		mat4 inverseViewTransform = inverse(viewTransform);
		vec4 transformedAxis = inverseViewTransform * vec4(0.0,1.0,0.0,0.0);

		mat4 newViewTransform = rotate(viewTransform, -0.5f*quarter_pi<float>(), vec3(transformedAxis));
		viewer()->setViewTransform(newViewTransform);
	}
	else if (key == GLFW_KEY_RIGHT && action == GLFW_RELEASE)
	{
		mat4 viewTransform = viewer()->viewTransform();
		mat4 inverseViewTransform = inverse(viewTransform);
		vec4 transformedAxis = inverseViewTransform * vec4(0.0, 1.0, 0.0, 0.0);

		mat4 newViewTransform = rotate(viewTransform, 0.5f*quarter_pi<float>(), vec3(transformedAxis));
		viewer()->setViewTransform(newViewTransform);
	}
	else if (key == GLFW_KEY_UP && action == GLFW_RELEASE)
	{
		mat4 viewTransform = viewer()->viewTransform();
		mat4 inverseViewTransform = inverse(viewTransform);
		vec4 transformedAxis = inverseViewTransform * vec4(1.0, 0.0, 0.0, 0.0);

		mat4 newViewTransform = rotate(viewTransform, -0.5f*quarter_pi<float>(), vec3(transformedAxis));
		viewer()->setViewTransform(newViewTransform);
	}
	else if (key == GLFW_KEY_DOWN && action == GLFW_RELEASE)
	{
		mat4 viewTransform = viewer()->viewTransform();
		mat4 inverseViewTransform = inverse(viewTransform);
		vec4 transformedAxis = inverseViewTransform * vec4(1.0, 0.0, 0.0, 0.0);

		mat4 newViewTransform = rotate(viewTransform, 0.5f*quarter_pi<float>(), vec3(transformedAxis));
		viewer()->setViewTransform(newViewTransform);
	}
	else if (key == GLFW_KEY_B && action == GLFW_RELEASE)
	{
		std::cout << "Starting benchmark" << std::endl;

		m_benchmark = true;
		m_startTime = glfwGetTime();
		m_frameCount = 0;
	}
	else if (key == GLFW_KEY_H && action == GLFW_RELEASE)
	{
		m_headlight = !m_headlight;
	}
	else if (key == GLFW_KEY_A && action == GLFW_RELEASE)
	{
		//add keyframe
		if (keyframes.size() >= 4)
		{
			std::cout << "error: maximum number of keyframes reached!\ndelete the last frame to create a new one!\n";
		}
		else
		{
			Keyframe kf;
			kf.backgroundColor = viewer()->backgroundColor();
			kf.explosionFactor = viewer()->m_explosion;
			kf.lightTransform = viewer()->lightTransform();
			kf.viewTransform = viewer()->viewTransform();
			keyframes.push_back(kf);
			std::cout << "Keyframe " << keyframes.size()<< " added\n";
		}

	}
	else if (key == GLFW_KEY_R && action == GLFW_RELEASE)
	{
		//remove keyframe
		if (keyframes.size() > 0)
		{
			keyframes.pop_back();
			std::cout << "Keyframe " << keyframes.size() + 1 << " removed\n";
		}
		else
		{
			std::cout << "error: no keyframes present\n";
		}
	}
	else if (key == GLFW_KEY_P && action == GLFW_RELEASE)
	{
		playAnimation();
	}
}

void CameraInteractor::playAnimation()
{
	//play animation
	if (keyframes.size() == 4)
	{
		std::vector<vec3> background_colors = { keyframes.at(0).backgroundColor, keyframes.at(1).backgroundColor, keyframes.at(2).backgroundColor, keyframes.at(3).backgroundColor };
		std::vector<float> explosion_factors = { keyframes.at(0).explosionFactor, keyframes.at(1).explosionFactor, keyframes.at(2).explosionFactor, keyframes.at(3).explosionFactor };
		std::vector<mat4> view_matrices = { keyframes.at(0).viewTransform, keyframes.at(1).viewTransform, keyframes.at(2).viewTransform, keyframes.at(3).viewTransform };
		std::vector<mat4> light_matrices = { keyframes.at(0).lightTransform, keyframes.at(1).lightTransform, keyframes.at(2).lightTransform, keyframes.at(3).lightTransform };

		unsigned long start_time = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
		unsigned long current_time = start_time;
		unsigned long end_time = start_time + 4000;
		while (current_time < end_time)
		{
			float idx = (current_time - start_time) / 4000.0f;
			viewer()->setBackgroundColor(cubicBezierVector(background_colors, idx));
			viewer()->m_explosion = cubicBezierScalar(explosion_factors, idx);
			viewer()->setViewTransform(cubicBezierMatrix(view_matrices, idx));
			viewer()->setLightTransform(cubicBezierMatrix(light_matrices, idx));
			viewer()->display();
			glfwSwapBuffers(viewer()->window());
			current_time = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
		}
	}
	else
	{
		std::cout << "error: not enough keyframes added\n";
	}
}

void CameraInteractor::mouseButtonEvent(int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		m_rotating = true;
		m_xPrevious = m_xCurrent;
		m_yPrevious = m_yCurrent;
	}
	else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
	{
		m_scaling = true;
		m_xPrevious = m_xCurrent;
		m_yPrevious = m_yCurrent;
	}
	else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS)
	{
		m_panning = true;
		m_xPrevious = m_xCurrent;
		m_yPrevious = m_yCurrent;
	}
	else
	{
		m_rotating = false;
		m_scaling = false;
		m_panning = false;
	}
}

void CameraInteractor::cursorPosEvent(double xpos, double ypos)
{
	m_xCurrent = xpos;
	m_yCurrent = ypos;

	if (m_light)
	{
		vec3 v = arcballVector(m_xCurrent, m_yCurrent);
		mat4 viewTransform = viewer()->viewTransform();

		mat4 lightTransform = inverse(viewTransform)*translate(mat4(1.0f), -0.5f*v*m_distance)*viewTransform;
		viewer()->setLightTransform(lightTransform);
	}

	if (m_rotating)
	{
		if (m_xCurrent != m_xPrevious || m_yCurrent != m_yPrevious)
		{
			vec3 va = arcballVector(m_xPrevious, m_yPrevious);
			vec3 vb = arcballVector(m_xCurrent, m_yCurrent);

			if (va != vb)
			{
				float angle = acos(max(-1.0f, min(1.0f, dot(va, vb))));
				vec3 axis = cross(va, vb);

				mat4 viewTransform = viewer()->viewTransform();
				mat4 lightTransform = viewer()->lightTransform();
				mat4 inverseViewTransform = inverse(viewTransform);
				vec4 transformedAxis = inverseViewTransform * vec4(axis, 0.0);

				mat4 newViewTransform = rotate(viewTransform, angle, vec3(transformedAxis));
				viewer()->setViewTransform(newViewTransform);

				if (m_headlight)
				{
					mat4 newLightTransform = rotate(lightTransform, angle, vec3(transformedAxis));
					viewer()->setLightTransform(newLightTransform);
				}
			}
		}

	}

	if (m_scaling)
	{
		if (m_xCurrent != m_xPrevious || m_yCurrent != m_yPrevious)
		{
			ivec2 viewportSize = viewer()->viewportSize();
			vec2 va = vec2(2.0f*float(m_xPrevious) / float(viewportSize.x) - 1.0f, -2.0f*float(m_yPrevious) / float(viewportSize.y) + 1.0f);
			vec2 vb = vec2(2.0f*float(m_xCurrent) / float(viewportSize.x) - 1.0f, -2.0f*float(m_yCurrent) / float(viewportSize.y) + 1.0f);
			vec2 d = vb - va;

			float l = std::abs(d.x) > std::abs(d.y) ? d.x : d.y;
			float s = 1.0;

			if (l > 0.0f)
			{
				s += std::min(0.5f, length(d));
			}
			else
			{
				s -= std::min(0.5f, length(d));
			}

			mat4 viewTransform = viewer()->viewTransform();
			mat4 newViewTransform = scale(viewTransform, vec3(s, s, s));
			viewer()->setViewTransform(newViewTransform);
		}
	}

	if (m_panning)
	{
		if (m_xCurrent != m_xPrevious || m_yCurrent != m_yPrevious)
		{
			ivec2 viewportSize = viewer()->viewportSize();
			float aspect = float(viewportSize.x) / float(viewportSize.y);
			vec2 va = vec2(2.0f*float(m_xPrevious) / float(viewportSize.x) - 1.0f, -2.0f*float(m_yPrevious) / float(viewportSize.y) + 1.0f);
			vec2 vb = vec2(2.0f*float(m_xCurrent) / float(viewportSize.x) - 1.0f, -2.0f*float(m_yCurrent) / float(viewportSize.y) + 1.0f);
			vec2 d = vb - va;

			mat4 viewTransform = viewer()->viewTransform();
			mat4 newViewTransform = translate(mat4(1.0), vec3(aspect*d.x, d.y, 0.0f))*viewTransform;
			viewer()->setViewTransform(newViewTransform);
		}
	}

	m_xPrevious = m_xCurrent;
	m_yPrevious = m_yCurrent;

}

void CameraInteractor::scrollEvent(double xoffset, double yoffset)
{
	mat4 viewTransform = viewer()->viewTransform();
	mat4 newViewTransform = translate(mat4(1.0), vec3(0.0f, 0.0f, (yoffset / 8.0)))*viewTransform;
	viewer()->setViewTransform(newViewTransform);
}

void CameraInteractor::display()
{
	if (m_benchmark)
	{
		m_frameCount++;

		mat4 viewTransform = viewer()->viewTransform();
		mat4 inverseViewTransform = inverse(viewTransform);
		vec4 transformedAxis = inverseViewTransform * vec4(0.0, 1.0, 0.0, 0.0);

		mat4 newViewTransform = rotate(viewTransform, pi<float>() / 180.0f, vec3(transformedAxis));
		viewer()->setViewTransform(newViewTransform);


		if (m_frameCount >= 360)
		{
			double currentTime = glfwGetTime();

			std::cout << "Benchmark finished." << std::endl;
			std::cout << "Rendered " << m_frameCount << " frames in " << (currentTime - m_startTime) << " seconds." << std::endl;
			std::cout << "Average frames/second: " << double(m_frameCount) / (currentTime - m_startTime) << std::endl;

			m_benchmark = false;
		}


	}

	if (ImGui::BeginMenu("Camera"))
	{
		static int projection = 0;
		ImGui::RadioButton("Perspective", &projection,0);
		ImGui::RadioButton("Orthographic", &projection,1);

		if ((projection == 0) != m_perspective)
		{
			m_perspective = (projection == 0);
			resetProjectionTransform();
		}

		ImGui::Checkbox("Headlight", &m_headlight);
		ImGui::EndMenu();
	}
}

void CameraInteractor::resetProjectionTransform()
{
	vec2 viewportSize = viewer()->viewportSize();
	framebufferSizeEvent(viewportSize.x, viewportSize.y);
}

void CameraInteractor::resetViewTransform()
{
	viewer()->setViewTransform(lookAt(vec3(0.0f, 0.0f, -m_distance), vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f)));
	viewer()->setLightTransform(lookAt(vec3(0.0f, 0.0f, -0.5f*m_distance), vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f)));
}

vec3 CameraInteractor::arcballVector(double x, double y)
{
	ivec2 viewportSize = viewer()->viewportSize();
	vec3 p = vec3(2.0f*float(x) / float(viewportSize.x)-1.0f, -2.0f*float(y) / float(viewportSize.y)+1.0f, 0.0);

	float length2 = p.x*p.x + p.y*p.y;

	if (length2 < 1.0f)
		p.z = sqrt(1.0f - length2);
	else
		p = normalize(p);

	return p;
}

float CameraInteractor::cubicBezierScalar(std::vector<float> control_points, float u)
{
	return pow((1 - u), 3) * control_points.at(0) + 3 * u * pow((1 - u), 2) * control_points.at(1) + 3 * pow(u, 2) * (1 - u) * control_points.at(2) + pow(u, 3) * control_points.at(3);
}

vec3 CameraInteractor::cubicBezierVector(std::vector<vec3> control_points, float u)
{
	return pow((1 - u), 3) * control_points.at(0) + 3 * u * pow((1 - u), 2) * control_points.at(1) + 3 * pow(u, 2) * (1 - u) * control_points.at(2) + pow(u, 3) * control_points.at(3);
}


mat4 CameraInteractor::cubicBezierMatrix(std::vector<mat4> control_points, float u)
{
	if (control_points.at(0) == control_points.at(1) && control_points.at(1) == control_points.at(2) && control_points.at(2) == control_points.at(3))
	{
		return control_points.at(0);
	}
	else
	{
		std::vector<vec3> translations;
		std::vector<vec3> scalars;
		std::vector<quat> rotations;
		for (int i = 0; i < control_points.size(); i++)
		{
			translations.push_back(extractTranslation(control_points.at(i)));
			scalars.push_back(extractScale(control_points.at(i)));
			rotations.push_back(quat_cast(control_points.at(i)));
		}
		vec3 translation = cubicBezierVector(translations, u);
		vec3 scale = cubicBezierVector(scalars, u);
		quat rotation = quaternionCubicSLERP(rotations, u);

		mat4 rotation_mat = mat4_cast(rotation);
		mat4 translation_mat = mat4(vec4(1.0f, 0.0f, 0.0f, 0.0f), vec4(0.0f, 1.0f, 0.0f, 0.0f), vec4(0.0f, 0.0f, 1.0f, 0.0f), vec4(translation, 1.0f));
		mat4 scaling_mat = mat4(0.0f);
		scaling_mat[0][0] = scale.x;
		scaling_mat[1][1] = scale.y;
		scaling_mat[2][2] = scale.z;
		scaling_mat[3][3] = 1.0f;

		return translation_mat * rotation_mat * scaling_mat;
	}
}

vec3 CameraInteractor::extractTranslation(mat4&  matrix)
{
	vec3 translation = matrix[3];
	matrix[3] = vec4(vec3(0.0f), 1.0f);
	return translation;
}

vec3 CameraInteractor::extractScale(mat4& matrix)
{
	vec3 x = matrix[0];
	vec3 y = matrix[1];
	vec3 z = matrix[2];
	float sx = length(x);
	float sy = length(y);
	float sz = length(z);

	matrix[0] = vec4((x * (1 / sx)), 0.0f);
	matrix[1] = vec4((y * (1 / sy)), 0.0f);
	matrix[2] = vec4((z * (1 / sz)), 0.0f);

	return vec3(sx, sy, sz);
}

quat CameraInteractor::quaternionCubicSLERP(std::vector<quat> quaternions, float u)
{
	quat q1 = mix(quaternions[0], quaternions[1], u);
	quat q2 = mix(quaternions[1], quaternions[2], u);
	quat q3 = mix(quaternions[2], quaternions[3], u);

	quat q4 = mix(q1, q1, u);
	quat q5 = mix(q2, q3, u);
	return mix(q4, q5, u);
}

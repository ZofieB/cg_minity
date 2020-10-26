#pragma once
#include "Interactor.h"
#include <glm/glm.hpp>
#include <vector>

namespace minity
{
	struct Keyframe
	{
		glm::vec3 backgroundColor;
		float explosionFactor = 0;
		glm::mat4 lightTransform;
		glm::mat4 viewTransform;
	};

	class Viewer;

	class CameraInteractor : public Interactor
	{
	public:
		CameraInteractor(Viewer * viewer);
		virtual void framebufferSizeEvent(int width, int height);
		virtual void keyEvent(int key, int scancode, int action, int mods);
		virtual void mouseButtonEvent(int button, int action, int mods);
		virtual void cursorPosEvent(double xpos, double ypos);
		virtual void scrollEvent(double xoffset, double yoffset);
		virtual void display();

		void resetProjectionTransform();
		void resetViewTransform();

		virtual float cubicBezierScalar(std::vector<float> control_points, float u);
		virtual glm::vec3 cubicBezierVector(std::vector<glm::vec3> control_points, float u);
		virtual glm::mat4 cubicBezierMatrix(std::vector<glm::mat4> control_points, float u);

	private:

		glm::vec3 arcballVector(double x, double y);
		glm::vec3 extractTranslation(glm::mat4& matrix);
		glm::vec3 extractScale(glm::mat4& matrix);
		glm::quat quaternionCubicSLERP(std::vector<glm::quat> quaternions, float u);

		void playAnimation();

		float m_fov = glm::radians(60.0f);
		float m_near = 0.125f;
		float m_far = 32768.0f;
		float m_distance = 2.0f*sqrt(3.0f);
		bool m_perspective = true;
		bool m_headlight = true;

		bool m_light = false;
		bool m_rotating = false;
		bool m_scaling = false;
		bool m_panning = false;
		bool m_benchmark = false;
		double m_startTime = 0.0;
		glm::uint m_frameCount = 0;
		double m_xPrevious = 0.0, m_yPrevious = 0.0;
		double m_xCurrent = 0.0, m_yCurrent = 0.0;

		std::vector<Keyframe> keyframes;
	};
}
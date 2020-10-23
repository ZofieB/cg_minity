#pragma once
#include "Renderer.h"
#include <memory>

#include <glm/glm.hpp>
#include <glbinding/gl/gl.h>
#include <glbinding/gl/enum.h>
#include <glbinding/gl/functions.h>

#include <globjects/VertexArray.h>
#include <globjects/VertexAttributeBinding.h>
#include <globjects/Buffer.h>
#include <globjects/Program.h>
#include <globjects/Shader.h>
#include <globjects/Framebuffer.h>
#include <globjects/Renderbuffer.h>
#include <globjects/Texture.h>
#include <globjects/base/File.h>
#include <globjects/TextureHandle.h>
#include <globjects/NamedString.h>
#include <globjects/base/StaticStringSource.h>

namespace minity
{
	class Viewer;

	class ModelRenderer : public Renderer
	{
	public:
		ModelRenderer(Viewer *viewer);
		virtual void display();

	private:

		std::unique_ptr<globjects::VertexArray> m_lightArray = std::make_unique<globjects::VertexArray>();
		std::unique_ptr<globjects::Buffer> m_lightVertices = std::make_unique<globjects::Buffer>();

		int m_shininess = 10;
		glm::vec3 m_ambient = glm::vec3(0.5, 0.5, 0.5);
		glm::vec3 m_specular = glm::vec3(1.0, 1.0, 1.0);
		glm::vec3 m_diffuse = glm::vec3(0.5, 0.5, 0.5);
		bool diff_tex = false;
		bool amb_tex = false;
		bool spec_tex = false;
		bool obj_norm_map = false;
		bool tan_norm_map = false;
		bool bump_map = false;
		bool sinus = false;
		bool sinus_sqr = false;
		float amplitude = 1.0f;
		float frequency = 200.0f;
		int m_explosion = 0;
	};

}
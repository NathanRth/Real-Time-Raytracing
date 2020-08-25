"use strict"

var phong_vert=`#version 300 es
uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;
uniform mat3 normalMatrix;

uniform highp sampler2D TUp;
uniform highp sampler2D TUn;
uniform highp usampler2D TUi;

out vec3 Po;
out vec3 No;

const uint TEXWIDTH = 4096u;

ivec2 indirect_indices(usampler2D indices, int is)
{
	uint i = uint(is);
	ivec2 i2 = ivec2(i%TEXWIDTH,i/TEXWIDTH);
	uint j = texelFetch(indices,i2,0).r;
	return ivec2(j%TEXWIDTH,j/TEXWIDTH);
}

void main()
{
	ivec2 ii = indirect_indices(TUi, gl_InstanceID*3 + gl_VertexID);
	vec3 Q = texelFetch(TUp,ii,0).rgb;
	vec3 N = texelFetch(TUn,ii,0).rgb;

	No = normalMatrix*N;
	vec4 P4 = viewMatrix * vec4(Q,1);
	Po = P4.xyz;
	gl_Position = projectionMatrix * P4;
}`;


var phong_store_frag=`#version 300 es
precision highp float;
in vec3 Po;
in vec3 No;
uniform vec3 light_pos;

layout(location=0) out vec3 position;
layout(location=1) out vec4 normal;

void main()
{
	position = Po;
	normal = vec4(No,2.0*float(gl_FrontFacing)-1.0); //w -1/1 back/front
}`;

var fs_vert = `#version 300 es
void main()
{
	uint id = uint(gl_VertexID);
	vec2 tc = vec2((id%2u), (id/2u)); // dans [0,1]
	vec2 p = vec2((id%2u), (id/2u))*2.0 - 1.0; // dans [-1,+1]
	gl_Position = vec4(p, 0, 1);
}
`;

var phong_def_frag=`#version 300 es
precision highp float;
uniform sampler2D TUp;
uniform sampler2D TUn;
uniform vec3 light_pos;

out vec4 frag_out;

void main()
{
	ivec2 coord = ivec2(gl_FragCoord.xy); // Attention ici texture de la taille de l'ecran
	vec3 P = texelFetch(TUp,coord,0).rgb;
	vec4 ND = texelFetch(TUn,coord,0).rgba;
	vec3 N = normalize(ND.xyz)*ND.w;
	vec3 L = normalize(light_pos-P);
	float lamb = 0.1+0.9*max(0.0,dot(N,L));
	frag_out = vec4(vec3(lamb),1);
}`;


var mesh_rend = Mesh.emptyRenderer();
var prg_phong_deferred_p1 = null;
var prg_phong_deferred_p2 = null;
var fbo_def = null;

var geom_tex_p = null;
var geom_tex_n = null;
var geom_tex_i = null;

var nb_tris = 0;
var nb_vert = 0;
var BB = null;
var mesh = null;

function loaded()
{
	BB = mesh.BB;

	if (geom_tex_p) { geom_tex_p.delete();}
	if (geom_tex_n) { geom_tex_n.delete();}
	if (geom_tex_i) { geom_tex_i.delete();}

	geom_tex_p = Texture2d();
	geom_tex_p.simple_params(gl.NEAREST);
	geom_tex_p.from_float_buffer(mesh.positions,3);

	geom_tex_n = Texture2d();
	geom_tex_p.simple_params(gl.NEAREST);
	geom_tex_n.from_float_buffer(mesh.normals,3);

	geom_tex_i = Texture2d();
	geom_tex_i.simple_params(gl.NEAREST);
	geom_tex_i.from_index_buffer(mesh.tris);

	nb_tris = mesh.tris.length / 3;
}

function resize_wgl(w,h)
{
	fbo_def.resize(w,h)
}

function init_wgl()
{
	prg_phong_deferred_p1 = ShaderProgram(phong_vert,phong_store_frag,'phong_p1');
	prg_phong_deferred_p2 = ShaderProgram(fs_vert,phong_def_frag,'phong_p2');
	mesh = Mesh.Wave(50);
	loaded();

	let texP = Texture2d();
	texP.simple_params(gl.NEAREST); // use texelfetch si no need to interpol.
	texP.init(gl.RGBA32F); // init only alloc in resize

	let texN = Texture2d();
	texN.simple_params(gl.NEAREST);
	texN.init(gl.RGBA32F);

	fbo_def = FBO_Depth([texP,texN]);

	ewgl.scene_camera.set_scene_radius(3);
	ewgl.scene_camera.set_scene_center(Vec3(0,0,0));

}

function draw_wgl()
{
	 // PASSE 1
	fbo_def.bind();
	gl.clearColor(0,0,0,1);
	gl.enable(gl.DEPTH_TEST);
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	const projection_matrix = ewgl.scene_camera.get_projection_matrix();
	const view_matrix = ewgl.scene_camera.get_view_matrix();

	prg_phong_deferred_p1.bind();
	Uniforms.viewMatrix =view_matrix;
	Uniforms.normalMatrix = view_matrix.inverse3transpose();
	Uniforms.projectionMatrix = projection_matrix;
	Uniforms.TUp = geom_tex_p.bind(0);
	Uniforms.TUn = geom_tex_n.bind(1);
	Uniforms.TUi = geom_tex_i.bind(2);
	gl.drawArraysInstanced(gl.TRIANGLES,0,3,nb_tris);

	unbind_fbo();

	// PASSE2
	prg_phong_deferred_p2.bind();
	Uniforms.TUp = fbo_def.texture(0).bind(0);
	Uniforms.TUn = fbo_def.texture(1).bind(1);
	gl.drawArrays(gl.TRIANGLE_STRIP,0,4);

	
}

ewgl.launch_3d();
 
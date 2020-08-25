"use strict"

var prg_rt_primary = null;

var random_tex = null;

var nb_samples = 3;
var nb_bounces = 5;
var ambiant_light = Vec3(0.1,0.1,0.1);
var ambiant_contribute = false;
var ambiant_intensity = 1.0;
var light_size = 0.5;
var light_intensity = 10;
var show_random = false;

function gen_random_tex()
{
    let size = window.innerWidth * window.innerHeight * 32;
    let random_vectors = create_Vec_buffer(3,size);
    for(let i = 0 ; i < size ; i++)
    {
        let x = Math.random();
        let y = Math.random();
        let phi = 2.0 * Math.PI * x;
        let theta = Math.acos(1.0 - 2.0*y);

        // uniform sampling on sphere
        var sinTheta = Math.sin(theta);
        var d = Vec3(Math.sin(phi)*sinTheta, Math.cos(phi)*sinTheta,Math.cos(theta));
        
        random_vectors.push(d);
    }
    random_tex = Texture3d();
    random_tex.alloc(window.innerWidth,window.innerHeight, 32 ,gl.RGB32F,random_vectors);
}

function init_wgl() {
    prg_rt_primary = ShaderProgramFromFiles('shaders/simple_raytrace.vert','shaders/simple_raytrace_v1.frag');
    ewgl.scene_camera.set_scene_radius(3);
    ewgl.scene_camera.set_scene_center(Vec3(0, 0, 0));
    ewgl.scene_camera.set_fov(50)
    console.warn("Computing random_vectors 3dTexture... it may take some time");
    gen_random_tex();

    UserInterface.begin("Interface", true, false);
    UserInterface.use_field_set('V',"(P) performance warning at high value");
    UserInterface.end_use();
    // render settings
    UserInterface.use_field_set('V',"Render settings");
    UserInterface.add_check_box("Continuous render (P)", ewgl.continuous_update, function(b){ewgl.continuous_update = b});
    UserInterface.add_slider("Samples (P)",1,40,nb_samples,function(e){nb_samples = e;},function(){return nb_samples;},0);
    UserInterface.add_slider("Max bounces",1,20,nb_bounces,function(e){nb_bounces = e;},function(){return nb_bounces;},0);
    
    UserInterface.end_use();

    UserInterface.use_field_set('V',"Ambiant light");
    UserInterface.add_check_box("Contribute", ambiant_contribute, function(b){ambiant_contribute = b});
    UserInterface.add_slider("Intensity",0,100,ambiant_intensity*10,function(e){ambiant_intensity = e/10.0;},function(){return ambiant_intensity;},2);
    UserInterface.add_slider("R",0,100,ambiant_light.x*100,function(e){ambiant_light.x = e/100.0;},function(){return ambiant_light.x;},2);
    UserInterface.add_slider("V",0,100,ambiant_light.y*100,function(e){ambiant_light.y = e/100.0;},function(){return ambiant_light.y;},2);
    UserInterface.add_slider("B",0,100,ambiant_light.z*100,function(e){ambiant_light.z = e/100.0;},function(){return ambiant_light.z;},2);
    UserInterface.end_use();

    UserInterface.use_field_set('V',"Light settings");
    UserInterface.add_slider("Size",0,100,light_size*100,function(e){light_size = e/100.0;},function(){return light_size;},2);
    UserInterface.add_slider("Intensity",0,100,light_intensity,function(e){light_intensity = e;},function(){return light_intensity*100;},0);
    UserInterface.end_use();

    UserInterface.use_field_set('V',"Debug");
    UserInterface.add_check_box("Show random values", show_random, function(b){show_random = b});
    UserInterface.end_use();

    UserInterface.end();

}

function draw_wgl() {
    let count = 0;
    gl.clearColor(0.1, 0, 0, 1);
    gl.enable(gl.DEPTH_TEST);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    const projection_matrix = ewgl.scene_camera.get_projection_matrix();
    const view_matrix = ewgl.scene_camera.get_view_matrix();

    prg_rt_primary.bind();

    Uniforms.NB_BOUNCES = nb_bounces;
    Uniforms.NB_SAMPLES = nb_samples;
    //Uniforms.NB_LIGHT_SAMPLES = nb_light_samples;

    Uniforms.random_map = random_tex.bind(0);
    Uniforms.show_random = show_random;

    Uniforms.ambiant_light = ambiant_light;
    Uniforms.ambiant_contribute = ambiant_contribute;
    Uniforms.ambiant_intensity = ambiant_intensity;

    //Uniforms.cam_pos = ewgl.scene_camera.s_center;
    Uniforms.viewMatrix = view_matrix;
    Uniforms.projectionMatrix = projection_matrix;
    Uniforms.winsize = [window.innerWidth, window.innerHeight];
    //Uniforms.winsize = [gl.canvas.clientWidth, gl.canvas.clientHeight];

    // cook torance             // red specular         // white diffuse        // light blue diffuse
    Uniforms.mat_albedo =       [Vec3(1.0, 0.0, 0.0),   Vec3(1.0, 1.0, 1.0),    Vec3(0.5, 0.5, 1.0)];
    Uniforms.mat_metallic =     [0.9,                   0.0,                    0.0];
    Uniforms.mat_roughness =    [0.4,                   1.0,                    1.0]; 
    
    Uniforms.NB_SPHERES = 3;
    Uniforms.spheres_pos =      [Vec3(2,1,-10), Vec3(5,1,-8),   Vec3(0,1,-5)    ];
    Uniforms.spheres_radius =   [1.0,           1.0,            1.0             ];
    Uniforms.spheres_matid =    [0,             1,              1               ];

    /*Uniforms.NB_PLANES = 1;
    Uniforms.planes_pos =       [Vec3(0,0,0),   Vec3(0,0,-13)];
    Uniforms.planes_normal =    [Vec3(0,1,0),   Vec3(0,0,1)];
    Uniforms.planes_matid =     [2,             2];*/

    Uniforms.NB_LIGHTS = 1;
    Uniforms.lights_pos =       [Vec3(3,3,-4),  Vec3(-4,3,-10)  ];
    Uniforms.lights_col =       [Vec3(1,1,1),   Vec3(1,1,1)     ];
    Uniforms.lights_intensity = [light_intensity*100,          10.0            ];
    Uniforms.lights_size =      [light_size,    0.5             ];
    
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    count++;
}

ewgl.launch_3d();
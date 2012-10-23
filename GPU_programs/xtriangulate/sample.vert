
uniform float rescale_normals;
uniform bool bumpy_surf;

varying vec3 N, v; 

void main()
{ 
    vec4 cc = gl_Vertex;
    vec3 nn = gl_Normal;
    
    if(bumpy_surf) {
        nn += vec3(sin(cc.x*30), sin(cc.y*30), sin(cc.z*30))*0.09;
        normalize(nn);
    }

    N = (gl_NormalMatrix * nn) * rescale_normals;
    v = vec3(gl_ModelViewMatrix * gl_Vertex);
        
    gl_FrontColor = gl_Color;
    gl_BackColor = gl_SecondaryColor;
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}


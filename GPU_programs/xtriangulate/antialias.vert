
uniform float rescale_normals;

//uniform vec3 light_pos[1];
varying vec3 N, v; 

void main()
{ 
    N = (gl_NormalMatrix * gl_Normal) * rescale_normals;
    v = vec3(gl_ModelViewMatrix * gl_Vertex);
    
    //L = normalize(gl_LightSource[0].position - V.xyz);
        
    gl_FrontColor = gl_Color;
    //gl_FrontSecondaryColor = gl_LightSource[0].ambient;

    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}


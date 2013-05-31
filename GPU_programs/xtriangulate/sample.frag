
varying vec3 N, v;

uniform bool do_depth_peeling;
uniform sampler2DRectShadow depth_tex;

void main() {

    if(do_depth_peeling) {
        // Lookup the previous layer's depth value at this fragment position
        vec4 depth = shadow2DRect(depth_tex, gl_FragCoord.stp);
        // If it is greater than or equal the current fragment's depth then
        // discard it.
        if(depth.x >= gl_FragCoord.z) { discard; }
    }
    
    vec4 clr = gl_Color;
    if(!gl_FrontFacing) {
        N = -N;
        clr.a = 0.55f;
    } else 
        clr.a = 0.55f;
  
    vec3 L = normalize(gl_LightSource[0].position.xyz - v); 
    vec3 E = normalize(-v); // we are in Eye Coordinates, so EyePos is (0,0,0)
    vec3 R = normalize(-reflect(L,N)); 

    float ndotl = max(dot(N, L), 0.0);
    vec3 color = (clr.rgb * ndotl) + gl_LightSource[0].ambient.rgb;

    float ndoth = max(dot(R, E), 0.0);
    ndoth = pow(ndoth, 32);
    color += vec3(ndoth);

    

    //vec4 col0 = texture2DRect(some_tex, gl_FragCoord.st);
    gl_FragColor = vec4(color, clr.a);
}

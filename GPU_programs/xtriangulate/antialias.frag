
uniform sampler2DRect framebuf;

uniform sampler2DRectShadow depth;

varying vec3 N, v;

void main() {
  
    /*vec4 sample[9];
    sample[0] = texture2DRect(tex, cc);
    sample[1] = texture2DRect(tex, gl_TexCoord[0].st + vec2(1, 0));
    sample[2] = texture2DRect(tex, gl_TexCoord[0].st + vec2(-1, 0));
    sample[3] = texture2DRect(tex, gl_TexCoord[0].st + vec2(0, 1));
    sample[4] = texture2DRect(tex, gl_TexCoord[0].st + vec2(0, -1));

    sample[5] = texture2DRect(tex, gl_TexCoord[0].st + vec2(-1, -1));
    sample[6] = texture2DRect(tex, gl_TexCoord[0].st + vec2(1, -1));
    sample[7] = texture2DRect(tex, gl_TexCoord[0].st + vec2(1, 1));
    sample[8] = texture2DRect(tex, gl_TexCoord[0].st + vec2(-1, 1));

//      gl_FragColor = (2.0*(sample[1] + sample[2] +  sample[3] + sample[4]) +
//          1.0*sample[0] + 6.0*(sample[5] +  sample[6] + sample[7] + sample[8])) / 16.0;*/

    vec3 L = normalize(gl_LightSource[0].position.xyz - v); 
    vec3 E = normalize(-v); // we are in Eye Coordinates, so EyePos is (0,0,0)
    vec3 R = normalize(-reflect(L,N)); 

    //vec3 NN = normalize(N);
    //vec3 NL = normalize(L);
    //vec3 NH = normalize(NL + vec3(0.0, 0.0, 1.0));

    float ndotl = max(dot(N, L), 0.0);
   
    vec3 color = (gl_Color.rgb * ndotl) + gl_LightSource[0].ambient.rgb;

    float ndoth = max(dot(R, E), 0.0);
    ndoth *= ndoth; // 2
    ndoth *= ndoth; // 4
    ndoth *= ndoth; // 8
    ndoth *= ndoth; // 16
    ndoth *= ndoth; // 32
   // ndoth *= ndoth; // 64
   // ndoth *= ndoth; // 128
    color += vec3(ndoth);
    
    //color = clamp(color, 0.0, 1.0);

    vec4 zbuf = shadow2DRect(depth, gl_FragCoord.stp);
    vec4 col = texture2DRect(framebuf, gl_FragCoord.st);

    vec4 color2;
   /* if(gl_FragCoord.z < zbuf.x) {
        // new fragment is in front of already drawn part
        color2 = mix(col, color, color.a);
        //discard;
    } else {
        // new fragment is behind
        color2 = mix(color, col, col.a);
    }*/

   
    color2 = vec4(color * col.rgb, 1.0);//mix(col, color, 0.9);

    //vec4 col0 = texture2DRect(some_tex, gl_FragCoord.st);
    gl_FragColor = vec4(color, gl_Color.a);
        //mix(vec4(1.0,0.0,0.0,1.0), gl_Color, 0.1); 
 
    // mix(a, b, t) = a * (1-t) + b * t
}


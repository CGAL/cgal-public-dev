uniform sampler2DRectShadow DT;

void main()
{ 
    // Lookup the previous layer's depth value at this fragment position
    vec4 depth = shadow2DRect( DT, gl_FragCoord.stp );

    // If it is greater than or equal the current fragment's depth then
    // discard it.
    if( depth.x >= gl_FragCoord.z ) { discard; }

    gl_FragColor = gl_Color;
}

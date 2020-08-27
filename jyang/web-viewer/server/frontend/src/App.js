import React, { Component } from 'react';
import * as THREE from "three";
import logo from './logo.svg';
import './App.css';

import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';
import io from 'socket.io-client';

let socket = io('http://127.0.0.1:3001');

// define scene
var scene = new THREE.Scene();
var camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
var renderer = new THREE.WebGLRenderer({ antialias: true });

var render = function () {
  renderer.render(scene, camera);
}

var onWindowResize = function () {
  console.log("RESIZE");
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
  render();
}

class App extends Component {
  componentDidMount() {
    this.ThreeJSInit();
    this.SocketIOInit();
    render();
  }

  ThreeJSInit() {
    // Init the camera
    camera.position.z = 2;

    // Init renderer
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.getElementById("widget").appendChild(renderer.domElement);

    // Init resizer
    window.addEventListener('resize', onWindowResize, false);

    // Init OrbitControls
    var controls = new OrbitControls(camera, renderer.domElement);
    controls.addEventListener('change', render); //use if there is no animation loop
    controls.minDistance = 0.1;
    controls.maxDistance = 10;
    controls.target.set(0, 0, 0);
    controls.update();
  }

  SocketIOInit() {
    socket.emit('message', 'Hello from React Frontend');

    socket.on('message', (message) => {
      console.log(message);
    });

    socket.on('vertices', (vertex_buffer) => {
      // add geometry
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertex_buffer, 3));
      var material = new THREE.PointsMaterial({ color: 0x00ff00, size: 0.003 });
      var points = new THREE.Points(geometry, material);
      scene.add(points);

      // render scene
      render();
    });

    socket.on('lines', (vertex_buffer) => {
      // add geometry
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertex_buffer, 3));
      var wireframe = new THREE.WireframeGeometry(geometry);
      var line = new THREE.LineSegments(wireframe);
      line.material.depthTest = false;
      line.material.opacity = 0.25;
      line.material.transparent = true;
      scene.add(line);

      // render scene
      render();
    });

    socket.on('triangles', (vertex_buffer) => {
      // add geometry
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertex_buffer, 3));
      var material = new THREE.MeshBasicMaterial({ color: 0x0000ff });
      var mesh = new THREE.Mesh(geometry, material);
      scene.add(mesh);

      // render scene
      render();
    });

    socket.on('default', (geometry_str) => {
      // add geometry
      var geometry = new THREE.Geometry();
      switch (geometry_str) {
        case 'cube': geometry = new THREE.BoxGeometry(); break;
        case 'sphere': geometry = new THREE.SphereGeometry(); break;
        default: geometry = new THREE.BoxGeometry(); break;
      }
      geometry.computeFaceNormals();
      var mesh = new THREE.Mesh(geometry, new THREE.MeshNormalMaterial());
      scene.add(mesh);

      // render scene
      render();
    });
  }

  render() {
    return (
      <div className="App">
        <div id="widget"/>
      </div>
    );
  }
}

export default App;

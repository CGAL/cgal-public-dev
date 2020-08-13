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
  constructor(props) {
    super(props);
  }

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
    controls.minDistance = 2;
    controls.maxDistance = 10;
    controls.target.set(0, 0, 0);
    controls.update();
  }

  SocketIOInit() {
    socket.emit('message', 'Hello from React Frontend');
    socket.on('message', (message) => {
      console.log(message);
    });
    socket.on('vertices', (vertices_str) => {
      // decode vertices
      var vertices = [];
      var xyz_str_list = vertices_str.trim().split(';');
      for (var xyz_str of xyz_str_list) {
        if (!xyz_str.length) continue;
        var xyz = xyz_str.trim().split(' ');
        vertices.push(parseFloat(xyz[0]), parseFloat(xyz[1]), parseFloat(xyz[2]));
      }

      // add geometry
      var geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
      var material = new THREE.PointsMaterial({color: 0x888888, size: 0.1});
      var points = new THREE.Points(geometry, material);
      scene.add(points)

      // render scene
      render();
    })
    socket.on('default', (geometry_str) => {
      var geometry = new THREE.Geometry();
      switch(geometry_str) {
        case 'cube': geometry = new THREE.BoxGeometry(); break;
        case 'sphere': geometry = new THREE.SphereGeometry(); break;
        default: geometry = new THREE.BoxGeometry(); break;
      }
      var material = new THREE.MeshBasicMaterial({ color: 0x00ff00 });
      var mesh = new THREE.Mesh(geometry, material);
      scene.add(mesh);
      render();
    })
  }

  render() {
    return (
      <div className="App">
        <div id="widget"></div>
        <header className="App-header">
          <img src={logo} className="App-logo" alt="logo" />
          <p>
            Edit <code>src/App.js</code> and save to reload.
          </p>
          <a
            className="App-link"
            href="https://reactjs.org"
            target="_blank"
            rel="noopener noreferrer"
          >
            Learn React
          </a>
        </header>
      </div>
    );
  }
}

export default App;

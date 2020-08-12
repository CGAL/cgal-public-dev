import React, { Component } from 'react';
import * as THREE from "three";
import logo from './logo.svg';
import './App.css';

import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls';
import io from 'socket.io-client';
let socket = io('http://127.0.0.1:3002');

class App extends Component {
  constructor(props) {
    super(props);
    this.state = {
      message: ''
    }
  }

  componentDidMount() {
    this.init();
    socket.on('message', (message) => {
      console.log(message);
    });
  }

  init() {
    // Create the scene
    var scene = new THREE.Scene();
    var camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    camera.position.z = 2;

    // Create renderer
    var renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.getElementById("widget").appendChild(renderer.domElement);

    // Create geometry
    var geometry = new THREE.BoxGeometry();
    var material = new THREE.MeshBasicMaterial( {color: 0x00ff00} );
    var cube = new THREE.Mesh(geometry, material);
    scene.add(cube);

    // Create render
    var render = function () {
      renderer.render(scene, camera);
    }
    render();

    // Create resizer
    var onWindowResize = function () {
      console.log("RESIZE");
      camera.aspect = window.innerWidth / window.innerHeight;
      camera.updateProjectionMatrix();
      renderer.setSize(window.innerWidth, window.innerHeight);
      render();
    }
    window.addEventListener('resize', onWindowResize, false);

    // Create OrbitControls
    var controls = new OrbitControls(camera, renderer.domElement);
    controls.addEventListener('change', render); //use if there is no animation loop
    controls.minDistance = 2;
    controls.maxDistance = 10;
    controls.target.set(0, 0, 0);
    controls.update();

  }

  render() {

    // const socket = socketIOClient('http://127.0.0.1:3002');
    // socket.on('message', (message) => {
    //   console.log(message);
    // });

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

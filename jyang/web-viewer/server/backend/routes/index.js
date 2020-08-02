var express = require('express');
const app = require('../app');
var router = express.Router();

// var http = require('http').createServer(app);
// var io = require('socket.io')(http);
var net = require('net');

/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'Express' });
});

// io.on('connection', (socket) => {
//   console.log('a user connected');
// })
net.createServer(function(connection) {
  console.log('client connected');
})

module.exports = router;

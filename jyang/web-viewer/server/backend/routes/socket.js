var net = require('net');

var server = net.createServer(function(connection) {
    console.log('client connected');
  });

server.listen(3001, '127.0.0.1');

module.exports = server;
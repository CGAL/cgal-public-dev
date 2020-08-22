var createError = require('http-errors');
var express = require('express');
var path = require('path');
var cookieParser = require('cookie-parser');
var logger = require('morgan');

var indexRouter = require('./routes/index');
var usersRouter = require('./routes/users');
// require('./routes/socket');

var app = express();

// view engine setup
app.set('views', path.join(__dirname, 'views'));
app.set('view engine', 'jade');

app.use(logger('dev'));
app.use(express.json());
app.use(express.urlencoded({ extended: false }));
app.use(cookieParser());
app.use(express.static(path.join(__dirname, 'public')));

app.use('/', indexRouter);
app.use('/users', usersRouter);

/* ------------------------------------------- */
//  communication between server backend (Express) and server fronend (React)
const server_react = require('http').createServer(app);
const io = require('socket.io')(server_react);

io.on('connection', (socket) => {
  io.emit('message', 'Hello from Express Backend');

  // connection info
  const host = socket.handshake.headers['host'].split(':');
  console.log('client connected from react at', host[0], ':', host[1]);

  // add 'message' event handler to this socket instance
  socket.on('message', (message) => {
    console.log(message);
  })
});

server_react.listen(3001, '127.0.0.1', () => {
  console.log('running react server on', server_react.address());
});

/* ------------------------------------------- */
// communication between cpp client and server backend
// Tcp packet size is 64K (65535 bytes), which is need to be considered
var server_cpp = require('net').createServer((socket) => {
    // connection info
    console.log('client connected from cpp at', socket.remoteAddress, ':', socket.remotePort);
    
    // add 'data' event handler to this socket instance
    socket.on('data', (data) => {
      console.log(socket.bytesRead, 'bytes', typeof data, 'data received from cpp:');
      console.log(data.buffer.byteLength)
      console.log(data.byteLength)
      
      // decode mode
      var mode; // first 4 bytes as mode
      switch (new DataView(data.buffer).getInt32(0)) {
        case 0: mode = 'vertices'; break;
        case 1: mode = 'lines'; break;
        case 2: mode = 'triangles'; break;
      }

      // decode elements
      var elements = [];
      for (var i = 4; i < 1060; i = i + 8) {
        elements.push(new DataView(data.buffer.slice(i, i+8)).getFloat64(0));
      }
      
      // resend the data to React Frontend
      console.log(mode, elements);
      io.emit(mode, elements);
    });

    var message = 'Hello from Express Backend';
    socket.end(message, () => {
      console.log(socket.bytesWritten, 'bytes', typeof message, 'data sent from Express backend:', message);
    });
  }).on('error', (err) => {
    // handle errors here.
    throw err;
  });

  server_cpp.listen(3002, '127.0.0.1', ()=> {
  console.log('running cpp server on', server_cpp.address());
});

/* ------------------------------------------- */
// catch 404 and forward to error handler
app.use(function(req, res, next) {
  next(createError(404));
});

// error handler
app.use(function(err, req, res, next) {
  // set locals, only providing error in development
  res.locals.message = err.message;
  res.locals.error = req.app.get('env') === 'development' ? err : {};

  // render the error page
  res.status(err.status || 500);
  res.render('error');
});

module.exports = app;

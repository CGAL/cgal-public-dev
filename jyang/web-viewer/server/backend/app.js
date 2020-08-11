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

// communication between client and server backend
var net = require('net');
var server = net.createServer(function(socket) {
    console.log('client connected from', socket.remoteAddress, ':', socket.remotePort);
    
    // set data encoding
    socket.setEncoding('utf-8');

    // add 'data' event handler to this socket instance
    socket.on('data', (data) => {
      console.log(socket.bytesRead, 'bytes', typeof data, 'data received:', data.toString('utf-8'));
    });

    var message = 'goodbye';
    socket.end(message, () => {
      console.log(socket.bytesWritten, 'bytes', typeof message, 'data sent:', message);
    });
  }).on('error', (err) => {
    // handle errors here.
    throw err;
  });

server.listen(3001, '127.0.0.1', ()=> {
  console.log('running server on', server.address());
});

//  try here
var http = require('http');
var socketIO = require('socket.io');
const server_react = http.createServer(app);
const io = socketIO(server_react);
io.on('connection', (socket) => {
  console.log('connected');
});

server_react.listen(3002, '127.0.0.1', () => {
  console.log('running server on', server_react.address());
});


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

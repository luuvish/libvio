#!/usr/bin/env node

var codec = require('./build/Release/codec');

codec.h264.main(process.argv[2], process.argv[3]);

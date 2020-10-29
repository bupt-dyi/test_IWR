/**
 * Local web server for the application.
 */

const opn			= require('opn');
const path		 	= require('path');
const args		    = require("optimist").argv;
const os			= require('os');
const fileutils     = require('../node_modules/ti_guicomposer_server/component-shares/server/util/fileutils.js');
const http 		 	= require('http');
const child_process = require('child_process');
const fs			= require('fs');
const serveStatic  	= require('serve-static');

/*******************************************************************
 * Get the app profile path
 *******************************************************************/
const getProfilePath = function(app) {
    const getWin     = function () {return path.join(os.homedir(), 'AppData', 'Local')};
    const getLinux   = function () {return path.join(os.homedir(), '.config')};
    const getOSX     = function () {return path.join(os.homedir(), 'Library', 'Application Support')};
    const getDefault = function () {return os.platform().startsWith('win') ? getWin() : getLinux()};

    let appDataPath = process.env['APPDATA'];
    if (!appDataPath) {
        switch (os.platform()) {
            case 'win32':  appDataPath = getWin(); break;
            case 'linux':  appDataPath = getLinux(); break;
            case 'darwin': appDataPath = getOSX(); break;
            default:       appDataPath = getDefault(); break;
        }
    }

    return path.join(appDataPath, app);
};

/* splash */
if (args.browser.indexOf('nw') > 0) {
	if (!/^darwin/.test(process.platform)) {
        const cleanup = function() {
            const profilePath = getProfilePath('splash' /* app/splash/package.json name property */);
            fileutils.rmdirSync(profilePath);
        };
        opn('', { app: [args.browser, path.join(__dirname, 'splash')], wait: false }).then(function(splashProcess) {
            splashProcess.on('exit', cleanup.bind(null,));
            splashProcess.on('SIGINT', cleanup.bind(null,));
            splashProcess.on('SIGUSR1', cleanup.bind(null,));
            splashProcess.on('SIGUSR2', cleanup.bind(null,));
            splashProcess.on('uncaughtException', cleanup.bind(null));
        });
    }
}

function AppServer() {
	this.keepAliveTimeout  	     = args.keepAliveTimeout || 120;
	this.port 					 = args.port || 0;
	this.keepAliveMonitorEnabled = args.keepAliveMonitorEnabled == null || args.keepAliveMonitorEnabled == 'true'
	this.browser				 = args.browser;
	this.baseDir				 = __dirname;
	this.cloudAgentDir			 = args.cloudAgentDir;
	this.cloudAgentPath			 = path.join(this.baseDir, "..", this.cloudAgentDir);

	process.chdir(this.baseDir);
}

/**
 * Start the CloudAgent host app.
 */
AppServer.prototype._startCloudAgentHostApp = function(callback) {
    const ext = (this.cloudAgentDir === 'win32' ? '.bat' : '.sh');
    const workingDir = path.join(this.cloudAgentPath, 'TICloudAgentHostApp');
    const fileToExec = path.join(workingDir, 'ticloudagent' + ext);

    this.cloudAgentHostApp = child_process.spawn(fileToExec, ['not_chrome'], { cwd: workingDir });
    this.cloudAgentHostApp.stdout.on('data', function (data) {
    	callback(null, data.toString());
	});
    this.cloudAgentHostApp.stderr.on('data', function (data) {
        callback(data.toString());
    });
}

/**
 * Stop the Server.
 */
AppServer.prototype.stop = function() {
    console.log('GUI Composer App server exit(0).')
    process.exit(0);
};

/**
 * Start the server.
 */
AppServer.prototype.start = function() {
	const self = this;

	const serve = serveStatic(self.baseDir, {index: ['index.html', 'index.htm']});
	const server = http.createServer(function(req, res) {
		function _endRequest (code, data, contentType) {
			res.writeHead(code, {'content-type': contentType || 'text/plain'});
			res.end(data);
		};

		/* On first request, start monitor the keep alive ping and shutdown the server if there is no ping */
		if (self.keepAliveMonitorEnabled && (self._keepAliveMonitorHandler == null)) {
			self._keepAlivePingTime = process.uptime();
			self._keepAliveMonitorHandler = setInterval(function() {
				if ((process.uptime() - self._keepAlivePingTime) >= self.keepAliveTimeout) {
					console.log('GUI Composer App Server has timeout!');
					self.stop();
				}
			}, 10000); // process last ping
		}

        /* public api */
        if (req.url.indexOf('/api') === 0) {
            var api = req.url.substring('/api'.length+1);
            switch (api) {
            case 'keepalive':
                self._keepAlivePingTime = process.uptime();
                console.log('GUI Composer App Server received keepalive ping request.');
                _endRequest(200, 'Server keep alive...');
                break;
            case 'shutdown':
                console.log('GUI Composer App Server received shutdown request.');
                _endRequest(200, 'Server shutting down...');
                if (args.browser.indexOf('nw') < 0) { /* only needs to shutdown for non-nodewebkit browser */
                    self.stop();
                }
                break;
            }

		/* runtime version.xml file */
		} else if (req.url.indexOf('version.xml') !=  -1) {
			var versionFile = fs.readFileSync(path.join(self.baseDir, '../version.xml'), 'utf-8');
			_endRequest(200, versionFile);

		/* server-config.js request */
		} else if (req.url.indexOf('server-config.js') != -1) {
			var config = 'window.gc = window.gc || {};\n';
			var configJson = {
				keepAliveMonitorEnabled: self.keepAliveMonitorEnabled,
				keepAliveTimeout: self.keepAliveTimeout,
				isOnline: false
			};
			config += 'window.gc.serverConfig = ' + JSON.stringify(configJson) + ';';
			_endRequest(200, config);

		/* agent.js route */
		} else if (req.url.indexOf('agent.js') != -1) {
			var agentJSPath = path.join(self.cloudAgentPath, "ticloudagent", "server", "public", "agent.js");
			var agentJS = fs.readFileSync(agentJSPath);
			_endRequest(200, agentJS);

		/* ticloudagent route */
		} else if (req.url.indexOf('/ticloudagent') === 0) {
			try {
				var reqPath = req.url.substring('/ticloudagent'.length+1);
				if (reqPath.indexOf('agent_config.js') === 0) {
					self._startCloudAgentHostApp(function(err, data) {
						if (err) {
							var agentConfig = {
								error: err
							}
							_endRequest(200, JSON.stringify(agentConfig));
						} else {
							var _data = JSON.parse(data);
							var agentConfig = {
								offline: true,
								version: _data.version,
								agentPort: _data.port
							};
							_endRequest(200, JSON.stringify(agentConfig));
						}
					});
				} else {
					var filePath = path.join(self.cloudagentDir, reqPath);
					_endRequest(200, fs.readFileSync(filePath, 'utf-8'))
				}

			} catch (e) {
				_endRequest(404, 'Bad Request: /ticloudagent');
			}

		/* scripts route */
		} else if (req.url.indexOf('/appscript') === 0) {
			if (args.script) {
				try {
					var filepath = args.script.trim();
					if (!path.isAbsolute(filepath)) {
						filepath = path.join(__dirname, filepath);
					}
					_endRequest(200, fs.readFileSync(filepath, 'utf-8'));

				} catch (e) {
					_endRequest(400, 'Failed to read app script file.');
				}
			} else {
				_endRequest(404, 'No app script defined');
			}

		/* log file route */
		} else if (req.url.indexOf('/logfile') === 0) {
			if (args.logfile) {
				var filepath = args.logfile.trim();
				if (!path.isAbsolute(filepath)) {
					filepath = path.join(__dirname, filepath);
				}
				_endRequest(200, filepath);

			} else {
				_endRequest(404, 'Log file not defined');
			}
		} else {
			serve(req, res, function() {
				_endRequest(404, 'File not found.');
			});
		}
	});

	server.listen(0, self.port, '127.0.0.1', function() {
		var port = server.address().port;
		console.log('GUI Composer App Server listening at port: ' + port);

		if (self.browser != null) {
			var url = 'http://127.0.0.1:' + port;

			if (self.browser.indexOf("nw") != -1) {
				/* app */
				/* workaround nw.exe can't accept dynamic port with package.json */ {
					var tmpDir = path.join(os.tmpdir(), '~gc-app-launcher');
					try {fs.mkdirSync(tmpDir)} catch (t) {};
					var packageJson = JSON.parse(fs.readFileSync(path.join(self.baseDir, 'package.json')));

					/* set the url include the port number */
					packageJson.main = url;

					/* copy the icon file - nw.exe require it to be relative to the package.json */
					if (packageJson.window) {
						var projectJson = JSON.parse(fs.readFileSync(path.join(self.baseDir, 'project.json')));
						if (projectJson.windowIcon) {
							var srcIconPath = path.join(self.baseDir, 'images', projectJson.windowIcon);
							var destIconPath = path.join(tmpDir, 'images', projectJson.windowIcon);

							if (!fs.existsSync(path.dirname(destIconPath)))
								fs.mkdirSync(path.dirname(destIconPath));

							fs.createReadStream(srcIconPath).pipe(fs.createWriteStream(destIconPath));
							packageJson.window.icon = 'images/' +  projectJson.windowIcon;
						}
					}

					/* write the package.json to the temp folder */
					fs.writeFileSync(path.join(tmpDir, 'package.json'), JSON.stringify(packageJson, undefined, 2));


                    /* clean up temp folder and terminated process */
                    const cleanup = function() {
                        fileutils.rmdirSync(tmpDir);
                        fileutils.removeItem(path.join(getProfilePath(packageJson.name), os.platform() === 'win32' ? 'User Data' : '', 'Default', 'Web Data'), self.stop);
                    };
                    process.on('exit', cleanup.bind(null,));
                    process.on('SIGINT', cleanup.bind(null,));
                    process.on('SIGUSR1', cleanup.bind(null,));
                    process.on('SIGUSR2', cleanup.bind(null,));
                    process.on('uncaughtException', cleanup.bind(null));

					/* launch nw */
					opn('', {app: [self.browser, tmpDir], wait: true }).then(cleanup.bind(self));
				}
			} else {

				opn(url, {app: self.browser, wait: true});
			}
		}
	});

};
new AppServer().start();
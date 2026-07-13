hexo.on('serverMiddleware', app => {
    app.use((req, res, next) => {
      res.setHeader('Content-Type', 'text/html; charset=utf-8');
      next();
    });
  });
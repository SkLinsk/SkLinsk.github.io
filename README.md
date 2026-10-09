# SkLinsk 学习笔记

一个直接托管于 GitHub Pages 的静态网站，不依赖 Hexo、Jekyll 或 npm 构建。

- `index.html`：笔记首页。
- `archives/index.html`：归档列表。
- `assets/site.css`：全站样式。
- `assets/theme.js`：白色、黑色与跟随系统模式。
- `notes/sqp/index.html`：SQP 阅读页；同目录保留 Markdown 与 PDF 下载。
- `notes/sqp/assets/katex/`：完整的公式样式、TeX 数学字体与许可证。
- `.nojekyll`：让 GitHub Pages 直接发布原始静态文件。

在 `master` 分支提交后，GitHub Pages 自动发布。添加笔记时，同时更新首页与归档中的入口。公式 HTML 已预先排版；修改公式需重新渲染对应 HTML。不要移除 KaTeX 的字体声明、定位和间距规则，也不要对数学 SVG 应用导航图标的全局样式。

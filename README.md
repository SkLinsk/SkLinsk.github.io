# SkLinsk 学习笔记

一个直接托管于 GitHub Pages 的静态网站，不依赖 Hexo、Jekyll 或 npm 构建。

- `index.html`：笔记首页。
- `archives/index.html`：归档列表。
- `assets/site.css`：全站样式。
- `assets/theme.js`：白色、黑色与跟随系统模式。
- `notes/sqp/index.html`：SQP 四篇系列总目录；同目录保留合订版 Markdown 与 PDF。
- `notes/sqp/basics/`、`curvature/`、`filter/`、`funnel/`：四篇独立阅读页，每篇都有 `note.md` 与 `note.pdf` 下载。
- `notes/sqp/assets/katex/`：完整的公式样式、TeX 数学字体与许可证。
- `notes/functional-analysis/`：泛函分析随笔《一个映射》，含 `note.md` 下载；公式复用 `notes/sqp/assets/katex/` 的样式与字体。
- `notes/cg/`：CG 与 Krylov 子空间学习笔记，含完整正交性推导、原数值例子、交互演示，以及 `note.md` / `note.pdf` 下载。
- `.nojekyll`：让 GitHub Pages 直接发布原始静态文件。

在 `master` 分支提交后，GitHub Pages 自动发布。添加笔记时，同时更新首页与归档中的入口。公式 HTML 已预先排版；修改公式需重新渲染对应 HTML。不要移除 KaTeX 的字体声明、定位和间距规则，也不要对数学 SVG 应用导航图标的全局样式。

## 网页写作后台

`/admin/` 提供 Markdown 写作、公式预览、图片、本机草稿和发布入口。正式 GitHub 登录需要部署 `auth/worker.mjs` 并配置自己的 GitHub App。按 `/admin/setup.html` 的步骤操作；密钥仅存于授权服务，不进入公开仓库。登录凭据只保存在当前页面内存。

文章由 `content/articles.json` 管理元数据，正文仍是各篇的 `note.md`。发布时同步生成 HTML、首页和归档，并检测并发改动。既有文章的交互组件会保留，旧 PDF 入口在正文改动后隐藏，避免内容版本不一致。本机草稿不会上传或跨设备同步。

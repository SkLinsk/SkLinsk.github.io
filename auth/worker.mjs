const json=(d,s=200)=>new Response(JSON.stringify(d),{status:s,headers:{'Content-Type':'application/json','Cache-Control':'no-store'}});
export default {async fetch(request,env){
 const url=new URL(request.url),origin=env.ALLOWED_ORIGIN||'https://sklinsk.github.io';
 if(url.pathname==='/health')return json({ready:!!(env.GITHUB_CLIENT_ID&&env.GITHUB_CLIENT_SECRET)});
 if(!env.GITHUB_CLIENT_ID||!env.GITHUB_CLIENT_SECRET)return json({error:'授权服务尚未配置 GitHub App'},503);
 if(url.pathname==='/auth'){
  const nonce=url.searchParams.get('nonce');if(!nonce||!/^[a-zA-Z0-9-]{20,100}$/.test(nonce))return json({error:'invalid nonce'},400);
  const state=crypto.randomUUID();const to=new URL('https://github.com/login/oauth/authorize');to.searchParams.set('client_id',env.GITHUB_CLIENT_ID);to.searchParams.set('redirect_uri',url.origin+'/callback');to.searchParams.set('state',state);to.searchParams.set('allow_signup','false');
  return new Response(null,{status:302,headers:{Location:to.href,'Set-Cookie':`writer_state=${state}.${nonce}; HttpOnly; Secure; SameSite=Lax; Path=/; Max-Age=600`,'Cache-Control':'no-store'}});
 }
 if(url.pathname!=='/callback')return json({error:'not found'},404);
 const saved=request.headers.get('Cookie')?.match(/(?:^|;\s*)writer_state=([^;]+)/)?.[1];const [expected,nonce]=saved?.split('.')||[];
 if(!expected||url.searchParams.get('state')!==expected)return json({error:'登录校验失效，请从写作页重新登录'},400);
 let payload;try{
  if(url.searchParams.get('error'))throw Error('GitHub 授权被取消');
  const code=url.searchParams.get('code');if(!code)throw Error('缺少授权码');
  const r=await fetch('https://github.com/login/oauth/access_token',{method:'POST',headers:{Accept:'application/json','Content-Type':'application/json'},body:JSON.stringify({client_id:env.GITHUB_CLIENT_ID,client_secret:env.GITHUB_CLIENT_SECRET,code,redirect_uri:url.origin+'/callback'})});const t=await r.json();if(!r.ok||!t.access_token)throw Error('GitHub 授权码交换失败');
  const headers={Accept:'application/vnd.github+json',Authorization:'Bearer '+t.access_token,'User-Agent':'sklinsk-writer-auth'};
  const u=await fetch('https://api.github.com/user',{headers}).then(r=>r.json());if(u.login?.toLowerCase()!==(env.ALLOWED_LOGIN||'SkLinsk').toLowerCase())throw Error('此后台只允许网站作者登录');
  const repo=await fetch('https://api.github.com/repos/SkLinsk/SkLinsk.github.io',{headers}).then(r=>r.json());if(!repo.permissions?.push)throw Error('请将 GitHub App 安装到网站仓库，并授予 Contents 读写权限');
  payload={type:'sklinsk-auth',nonce,token:t.access_token};
 }catch(e){payload={type:'sklinsk-auth',nonce,error:e.message}}
 const data=JSON.stringify(payload).replace(/</g,'\\u003c'),target=JSON.stringify(origin);const html=`<!doctype html><meta charset="utf-8"><meta name="referrer" content="no-referrer"><title>GitHub 登录</title><p>登录处理完成，可以关闭此窗口。</p><script>if(window.opener){window.opener.postMessage(${data},${target});window.close()}</script>`;
 return new Response(html,{headers:{'Content-Type':'text/html; charset=utf-8','Cache-Control':'no-store','Referrer-Policy':'no-referrer','Content-Security-Policy':"default-src 'none'; script-src 'unsafe-inline'; frame-ancestors 'none'; base-uri 'none'",'Set-Cookie':'writer_state=; HttpOnly; Secure; SameSite=Lax; Path=/; Max-Age=0'}});
}};

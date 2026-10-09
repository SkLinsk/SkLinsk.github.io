(()=>{
 const key='sklinsk-color-mode',root=document.documentElement;
 const system=matchMedia('(prefers-color-scheme: dark)');
 let mode='system';try{mode=localStorage.getItem(key)||'system'}catch{}
 if(!['light','dark','system'].includes(mode))mode='system';
 function apply(){root.dataset.theme=mode==='system'?(system.matches?'dark':'light'):mode;document.querySelectorAll('[data-mode]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.mode===mode)))}
 apply();
 document.addEventListener('DOMContentLoaded',()=>{
  apply();document.querySelectorAll('[data-mode]').forEach(b=>b.addEventListener('click',()=>{mode=b.dataset.mode;try{localStorage.setItem(key,mode)}catch{}apply()}));
  const toc=document.querySelector('.toc');if(toc){const media=matchMedia('(max-width:900px)');const sync=()=>{toc.open=!media.matches};sync();media.addEventListener('change',sync)}
 });
 system.addEventListener('change',apply);
 addEventListener('storage',e=>{if(e.key===key){mode=['light','dark','system'].includes(e.newValue)?e.newValue:'system';apply()}});
})();

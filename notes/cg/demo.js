(()=>{
 const panel=document.querySelector('.cg-demo');if(!panel)return;
 const steps=[
  {u:'(0, 0)',r:'(1, 2)',norm:'√5 ≈ 2.236068',note:'初始点：p₀ = r₀，开始第一次精确线搜索。'},
  {u:'(1/4, 1/2)',r:'(−1/2, 1/4)',norm:'√5/4 ≈ 0.559017',note:'α₀ = 1/4；新残差与初始残差正交。'},
  {u:'(1/11, 7/11)',r:'(0, 0)',norm:'0',note:'α₁ = 4/11；到达精确解，残差为零。'}
 ];
 let k=0;const prev=panel.querySelector('[data-cg-prev]'),next=panel.querySelector('[data-cg-next]');
 function render(){for(let j=1;j<=2;j++){panel.querySelector('#cg-leg-'+j).style.visibility=j<=k?'visible':'hidden';panel.querySelector('#cg-point-'+j).style.visibility=j<=k?'visible':'hidden';if(j===1)panel.querySelector('#cg-label-'+j).style.visibility=j<=k?'visible':'hidden';}panel.querySelector('[data-cg-k]').textContent=k;panel.querySelector('[data-cg-u]').textContent=steps[k].u;panel.querySelector('[data-cg-r]').textContent=steps[k].r;panel.querySelector('[data-cg-norm]').textContent=steps[k].norm;panel.querySelector('[data-cg-note]').textContent=steps[k].note;prev.disabled=k===0;next.disabled=k===2;}
 prev.addEventListener('click',()=>{k=Math.max(0,k-1);render()});next.addEventListener('click',()=>{k=Math.min(2,k+1);render()});panel.querySelector('[data-cg-reset]').addEventListener('click',()=>{k=0;render()});render();
})();

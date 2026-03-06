// Hide TOC for specific articles where tabsets break TOC rendering
(function() {
  var nameEl = document.querySelector('.template-article .name code');
  if (nameEl) {
    var noTocArticles = ['pbmc_vignette.Rmd'];
    if (noTocArticles.indexOf(nameEl.textContent.trim()) !== -1) {
      document.body.classList.add('no-toc');
    }
  }
})();

const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const p = await b.newPage();
  await p.emulateMedia({ colorScheme: 'light', media: 'print' });
  await p.goto('file://' + process.argv[2], { waitUntil: 'networkidle', timeout: 60000 });
  await p.evaluate(() => document.fonts.ready); await p.waitForTimeout(800);
  console.log('font check:', await p.evaluate(() => document.fonts.check('11pt "STIX Two Text"')));
  await p.pdf({ path: process.argv[3], format: 'Letter', printBackground: true, preferCSSPageSize: true,
    displayHeaderFooter: true, headerTemplate: '<span></span>',
    footerTemplate: '<div style="font-size:8pt;width:100%;text-align:center;font-family:serif;color:#333"><span class="pageNumber"></span></div>' });
  await b.close();
})();

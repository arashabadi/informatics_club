import path from "node:path";

const ROOT = "/Users/ashi/github/informatics_club";
const ARTIFACT_TOOL =
  "/Users/ashi/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";
const SOURCE = path.join(
  ROOT,
  "sctalk/06_seurat-vs-scanpy/06_sctalk0906_ready_to_read_notes.pptx",
);
const OUT = path.join(
  ROOT,
  "sctalk/06_seurat-vs-scanpy/06_sctalk0906_ready_to_read_notes_farsi.pptx",
);

const { PresentationFile, FileBlob } = await import(ARTIFACT_TOOL);

const FARSI = {
  7: `این اسلاید بخش اول Figure 1 را توضیح میدهد و نقطه شروع اختلاف Seurat و Scanpy را نشان میدهد. در پنل A سه نوع شیء مقایسه شده است: سلولها، ژنها، و ژنهای بسیار متغیر یا HVG. نکته مهم این است که اختلاف اصلی از همان انتخاب HVG شروع میشود، نه از این که Scanpy اصلا ژن نداشته باشد. عدد صفر در Aii به معنی «ژن صفر در Scanpy» نیست؛ یعنی در آن مقایسه، ژن اختصاصی Scanpy خارج از مجموعه Seurat وجود ندارد، چون ژنهایی که Scanpy برای ادامه نگه داشته داخل مجموعه بزرگتر Seurat قرار گرفته اند. بنابراین داستان این پنل این است که دو پکیج از نظر ورودی نهایی برای تحلیل پایین دست، مخصوصا HVGها، مسیر متفاوتی میسازند.

برای ارائه، این اسلاید را باید به عنوان نقطه آغاز زنجیره اختلافها توضیح داد. وقتی HVGها فرق کنند، PCA، گراف همسایگی، خوشه بندی، UMAP و marker geneها هم میتوانند عوض شوند. پس پیام شکل این نیست که یکی از ابزارها اشتباه ساده ای کرده؛ پیام این است که defaultهای به ظاهر معمول در ابتدای pipeline میتوانند مسیر کل تحلیل single-cell را تغییر بدهند.`,
  8: `در این بخش Figure 1، اختلاف از انتخاب ژن وارد فضای PCA و سپس گراف همسایگی میشود. پنل B نشان میدهد که شکل کلی PCA هنوز قابل مقایسه است، اما مقدار واریانس توضیح داده شده و جهت eigenvectorها یکسان نیست. این یعنی دو workflow داده را کاملا به یک مختصات کم بعدی تبدیل نکرده اند. در پنل C همین تفاوت وارد KNN/SNN graph میشود؛ نمودارهای Jaccard و degree ratio نشان میدهند که بسیاری از سلولها همسایه های مشترک مشابهی در دو پکیج ندارند.

برای کسی که مقاله را نخوانده، نکته کلیدی این است که clustering فقط روی expression خام انجام نمیشود؛ روی گرافی انجام میشود که بعد از PCA ساخته شده است. اگر گراف عوض شود، خوشه ها و UMAP هم ممکن است عوض شوند. پس این اسلاید مرحله میانی زنجیره را نشان میدهد: اختلافهای کوچک در feature selection حالا تبدیل به تفاوتهای ساختاری در فضای PCA و شبکه سلولی شده اند.`,
  9: `این اسلاید نتیجه قابل دیدن همان اختلافهای upstream را در پنلهای D و E نشان میدهد. در پنل D دو UMAP برای Seurat و Scanpy دیده میشود و نمودار alluvial وسط نشان میدهد که cluster assignmentها چقدر بین دو روش جابه جا شده اند. ساختار کلی شاید شبیه باشد، اما نگاشت clusterها کاملا یک به یک نیست. پنل E وارد differential expression میشود و overlap ژنهای marker معنی دار را بین دو workflow مقایسه میکند.

داستان ارائه این است که اختلاف نرم افزاری فقط در یک عدد فنی پنهان نمیماند. در نهایت کسی که تحلیل را میخواند، UMAP، clusterها و marker geneها را میبیند و بر اساس آنها تفسیر زیستی میکند. این شکل نشان میدهد defaultهای Seurat و Scanpy میتوانند همان خروجیهایی را تغییر دهند که معمولا در مقاله یا جلسه journal club درباره شان بحث میشود.`,
  10: `Figure 2 یک مقیاس برای فهمیدن بزرگی اختلاف Seurat و Scanpy میسازد. نویسنده ها پرسیده اند اگر فقط داخل یک پکیج داده را کم کنیم، مثلا readها یا cellها را downsample کنیم، چه مقدار کاهش داده لازم است تا تفاوت به اندازه تفاوت Seurat در برابر Scanpy برسد. در پنل A محور افقی fraction خوانشها را نشان میدهد و در پنل B fraction سلولها را. رنگ نارنجی Seurat و آبی Scanpy است و خط مرجع تفاوت بین دو پکیج روی داده کامل را نشان میدهد.

پیام مهم برای ارائه این است که انتخاب پکیج میتواند از نظر اثرگذاری هم اندازه یک perturbation بزرگ در داده باشد. در read downsampling، بسیاری از معیارها وقتی به کمتر از حدود پنج درصد read میرسند به سطح اختلاف بین پکیجها نزدیک میشوند. در cell downsampling این آستانه معمولا بالاتر است، ولی هنوز نشان میدهد تفاوت defaultها از نظر عملی کوچک نیست.`,
  11: `Figure 3 نشان میدهد مسئله reproducibility فقط Seurat در برابر Scanpy نیست؛ نسخه های نرم افزار هم مهم هستند. ردیف بالا Seurat v5 را با v4 مقایسه میکند، ردیف وسط Scanpy v1.9 را با v1.4، و ردیف پایین Cell Ranger v7 را با v6. پنلها به ترتیب تعداد markerهای معنی دار، log fold-change و adjusted p-value را بررسی میکنند. چون همه اینها مستقیم روی differential expression اثر دارند، خروجی marker geneها میتواند با تغییر نسخه عوض شود.

برای ارائه باید تاکید کرد که این اسلاید یک هشدار روش شناسی است. وقتی در مقاله نوشته میشود «از Seurat استفاده کردیم» یا «از Scanpy استفاده کردیم»، این کافی نیست. نسخه نرم افزار، تنظیمات محاسبه logFC، روش p-value adjustment، و حتی نسخه Cell Ranger قبل از ورود به تحلیل single-cell میتوانند نتیجه نهایی DE را تغییر دهند.`,
  12: `Figure 4 اختلافهای computational را به پیام زیستی وصل میکند. در پنل A خروجی CellTypist نشان میدهد که بخشی از سلولها بین pipeline پیش فرض Seurat و Scanpy annotation متفاوت میگیرند، مخصوصا در زیرگروههای monocyte و T cell. پنلهای B و C برای cluster 15 volcano plotهای Seurat و Scanpy را نشان میدهند؛ رنگ قرمز ژنهای upregulated معنی دار، آبی ژنهای downregulated معنی دار، و خاکستری ژنهای غیر معنی دار است. پنل D overlap ژنهای biologically significant را خلاصه میکند و پنل E نتیجه pathway analysis با Enrichr را نشان میدهد.

این اسلاید برای کسی که مقاله را نخوانده، بخش «چرا مهم است؟» را جواب میدهد. اگر فقط clusterها یا marker geneها کمی فرق کنند، ممکن است فکر کنیم مسئله فنی است. اما وقتی cell type assignment، ژنهای DE و pathwayها هم تغییر میکنند، تفسیر زیستی dataset عوض میشود. یعنی defaultهای نرم افزاری میتوانند به نتیجه علمی متفاوت برسند.`,
  13: `این اسلاید نقشه کدهای مقاله را نشان میدهد. سمت matrix_generation مربوط به مسیر پایین دستی تر از FASTQ تا count matrix است: configها خوانده میشوند، داده ها دانلود یا downsample میشوند، و ماتریس با kb یا Cell Ranger ساخته میشود. سمت analysis لایه ای است که شکلهای مقاله را تولید میکند؛ فایلهای YAML مشخص میکنند هر شکل با چه ورودی و چه تنظیماتی ساخته شود، و R Markdownها و scriptهای کمکی plotting، statistics، differential expression و رفتارهای نسخه های مختلف Scanpy را مدیریت میکنند.

برای آماده شدن ارائه، باید این repo را به دو سطح ذهنی تقسیم کرد. اگر فقط میخواهید بفهمید یک شکل چگونه ساخته شده، اول سراغ YAML همان شکل و Rmd مربوطه بروید. اگر میخواهید count matrixها را از raw FASTQ بازسازی کنید، آن وقت وارد matrix_generation میشوید. این جداسازی نشان میدهد مقاله سعی کرده هم upstream matrix generation و هم downstream analysis را قابل ردیابی کند.`,
  17: `Supplemental Figure 1 نقشه کلی مطالعه است. ردیف بالا منابع مختلف variability را نشان میدهد: عمق sequencing، تعداد سلول، نرم افزار ساخت ماتریس، پکیج تحلیل، نسخه نرم افزار، و stochasticity یا اثر random seed. ردیف پایین pipeline تحلیل single-cell را از filtering و normalization تا HVG selection، PCA، ساخت گراف، clustering، UMAP و differential expression دنبال میکند. رنگ قرمز منبع تغییرات را مشخص میکند و رنگ آبی metricهایی را که برای مقایسه خروجیها استفاده شده اند.

این اسلاید برای شروع بخش supplemental خیلی مهم است چون نشان میدهد مقاله فقط یک مقایسه موردی انجام نداده است. نویسنده ها میخواهند بفهمند در کدام مرحله های workflow، نتیجه میتواند از مسیر اصلی منحرف شود. بنابراین این شکل را میتوان به عنوان نقشه راه خواند: هر supplemental figure بعدی یکی از همین منبع های variability یا یکی از همین checkpointهای pipeline را باز میکند.`,
  18: `Supplemental Figure 2 میپرسد اگر به Scanpy آرگومانهای شبیه Seurat بدهیم، آیا خروجیها شبیه تر میشوند یا نه. اینجا تحلیل sequential است، یعنی خروجی هر مرحله وارد مرحله بعد میشود و اگر در اول pipeline اختلافی ایجاد شود، همان اختلاف به PCA، گراف، clustering، UMAP و DE منتقل میشود. پنل A filtering و HVG را مقایسه میکند، پنل B PCA را، پنل C گراف KNN/SNN را، پنل D cluster و UMAP را، و پنل E marker geneها را.

نکته مهم این است که تنظیم کردن argumentها کمک میکند، اما مثل دکمه reset عمل نمیکند. اگر early objectها مثل HVGها متفاوت باشند، حتی تنظیمات بعدی هم روی داده دقیقا یکسان اجرا نمیشوند. پس این شکل نشان میدهد harmonization واقعی بین Seurat و Scanpy سخت تر از این است که فقط چند parameter را شبیه کنیم.`,
  19: `Supplemental Figure 3 طراحی controlled analysis را نشان میدهد. در این طراحی، بعد از هر مرحله، object ساخته شده با Scanpy با object ساخته شده توسط Seurat جایگزین میشود تا مرحله بعدی در هر دو مسیر ورودی یکسان ببیند. بنابراین اگر در یک پنل اختلاف دیده شود، احتمال بیشتری دارد که آن اختلاف از خود function یا default همان مرحله آمده باشد، نه از خطا یا اختلاف انباشته شده در مراحل قبلی.

برای ارائه، این شکل را میتوان به عنوان آزمایش تشخیصی توضیح داد. Figure 1 میگوید اختلافها در workflow واقعی انباشته میشوند، اما این supplement میپرسد اگر accumulation را کنترل کنیم، کدام مرحله ها هنوز متفاوت میمانند. پنلهای filtering، PCA، گراف، clustering، UMAP و DE به ترتیب نشان میدهند کدام بخشها از اختلاف ناشی از local implementation هستند.`,
  20: `Supplemental Figure 4 همان controlled analysis را با Scanpy دارای آرگومانهای Seurat-like ترکیب میکند. یعنی از یک طرف ورودی هر مرحله کنترل شده است، و از طرف دیگر obvious parameterها تا حد امکان هماهنگ شده اند. پنلها دوباره همان مسیر آشنا را از filtering و HVG تا PCA، SNN graph، clustering، UMAP و differential expression نشان میدهند.

این شکل برای فهمیدن مرز harmonization مهم است. اگر در جایی agreement بهتر شود، یعنی parameter choice نقش اصلی داشته است. اگر با وجود ورودی کنترل شده و آرگومانهای aligned هنوز اختلاف بماند، یعنی مسئله عمیق تر است: ممکن است implementation، algorithm، یا default بدون معادل دقیق در دو پکیج متفاوت باشد.`,
  21: `Supplemental Figure 5 جهت هماهنگ سازی را برعکس میکند: این بار Seurat تا حد امکان Scanpy-like تنظیم میشود، در حالی که ورودی هر مرحله هم کنترل شده است. این تقارن مهم است چون اجازه نمیدهد یکی از پکیجها را «حقیقت پایه» فرض کنیم. پنلها همان مرحله های filtering، PCA، گراف، clustering، UMAP و DE را در شرایط کنترل شده بررسی میکنند.

پیام این اسلاید این است که reproducibility بین ابزارها فقط ترجمه یک لیست parameter از یک زبان به زبان دیگر نیست. بعضی تنظیمات معادل کامل ندارند، بعضی objectها در دو پکیج متفاوت ساخته یا ذخیره میشوند، و جهت harmonization خودش میتواند نتیجه را عوض کند. پس اگر بخواهیم تحلیل را بین Seurat و Scanpy بازتولید کنیم، باید architecture هر دو workflow را بفهمیم.`,
  22: `Supplemental Figure 6 روی UMAP-derived graph تمرکز میکند و نشان میدهد UMAP اینجا فقط یک تصویر تزئینی در آخر pipeline نیست. هر ردیف یک حالت مقایسه است: default یا aligned، و sequential یا controlled. در هر ردیف، density plotها overlap همسایگی را نشان میدهند، alluvial plotها Leiden clusterها را مقایسه میکنند، و UMAPها ساختار clusterها را روی embedding نشان میدهند.

این اسلاید به presenter کمک میکند توضیح دهد چرا اختلاف UMAP مهم است. اگر UMAP-derived neighborhood یا clusterها تغییر کنند، تفاوت فقط در شکل ظاهری نقاط نیست؛ ساختار محلی سلولها هم متفاوت خوانده شده است. بنابراین UMAP در این مقاله یک خروجی نمایشی ساده نیست، بلکه بخشی از مقایسه ساختاری بین workflowهاست.`,
  23: `Supplemental Figure 7 Seurat را با خودش مقایسه میکند وقتی readها به چهار درصد داده اصلی کاهش داده شده اند. چون هر دو طرف مقایسه Seurat هستند، اختلافها از انتخاب پکیج نمی آیند بلکه از کاهش read depth می آیند. پنل A overlap سلولها، ژنها و HVGها را نشان میدهد؛ پنلهای B و C اثر را در PCA و گراف همسایگی دنبال میکنند؛ پنل D cluster و UMAP را مقایسه میکند؛ و پنل E marker geneها را بررسی میکند.

این شکل یک calibration برای Figure 2 است. یعنی قبل از این که بگوییم تفاوت Seurat و Scanpy بزرگ است، باید ببینیم یک perturbation شدید داده داخل یک پکیج چه شکلی دارد. چهار درصد read downsampling یک تغییر شدید است، و این supplement نشان میدهد چنین تغییر شدیدی در checkpointهای مختلف pipeline چه اثری میگذارد.`,
  24: `Supplemental Figure 8 همان آزمایش چهار درصد read downsampling را داخل Scanpy انجام میدهد. ساختار پنلها شبیه Seurat است: filtering و HVG، سپس PCA، گراف، clustering و UMAP، و در آخر differential expression. چون package ثابت مانده، تغییرات نشان دهنده اثر از دست دادن readها در workflow Scanpy است، نه اختلاف بین دو ابزار.

برای ارائه، این اسلاید را باید کنار Supplementary Figure 7 دید. سوال این است که وقتی فقط read depth کم میشود، Seurat و Scanpy هر کدام چقدر self-consistent میمانند. این مقایسه کمک میکند اثر data depth را از اثر software default جدا کنیم و بفهمیم اختلاف اصلی مقاله در چه مقیاسی قرار میگیرد.`,
  25: `Supplemental Figure 9 به جای readها، تعداد سلولها را کاهش میدهد. در اینجا Seurat با خودش مقایسه میشود، اما فقط شانزده درصد سلولها نگه داشته شده اند. پنل A نشان میدهد مجموعه سلولها، ژنها و HVGها چطور تغییر میکنند. پنلهای B و C اثر را در PCA و graph neighborhood دنبال میکنند، پنل D cluster و UMAP را نشان میدهد، و پنل E marker geneها را بررسی میکند.

نکته این اسلاید این است که کم کردن سلولها نوع دیگری از perturbation است. در read downsampling هر سلول اطلاعات کمتری دارد، اما در cell downsampling ترکیب نمونه و representation جمعیت سلولی تغییر میکند. بنابراین این شکل نشان میدهد کم شدن تعداد سلولها چگونه میتواند ساختار خوشه ها و markerها را حتی داخل یک workflow ثابت تغییر دهد.`,
  26: `Supplemental Figure 10 نسخه Scanpy همان cell downsampling شانزده درصدی است. دوباره package ثابت است و فقط تعداد سلولها کم شده، پس تغییرها به sample-size loss مربوط هستند. پنلها از filtering و HVG شروع میکنند، بعد وارد PCA و graph construction میشوند، سپس clustering و UMAP را نشان میدهند، و در پایان differential expression را مقایسه میکنند.

برای کسی که مقاله را آماده ارائه میکند، این اسلاید همراه Supplementary Figure 9 یک جفت مقایسه ای است. باید دید وقتی تعداد سلولها شدیدا کم میشود، آیا Scanpy مثل Seurat رفتار میکند یا حساسیت متفاوتی دارد. این به نویسنده ها کمک میکند تفاوت پکیجها را در کنار تفاوت ناشی از اندازه نمونه تفسیر کنند.`,
  27: `Supplemental Figure 11 خلاصه کمی read downsampling را در چندین metric نشان میدهد. هر پنل یک جنبه pipeline را دنبال میکند: retained cells، retained genes، HVG overlap، PCA، graph degree، clustering agreement، UMAP-neighborhood overlap، marker gene overlap، log fold-change و عبور از threshold adjusted p-value. رنگ نارنجی Seurat و آبی Scanpy است و خط dashed تفاوت Seurat در برابر Scanpy روی داده کامل را نشان میدهد.

طرز خواندن این شکل این است که خط dashed را benchmark در نظر بگیریم. هر جا curve downsampling به آن خط برسد، اثر کاهش read داخل یک پکیج به اندازه اثر تغییر package شده است. بنابراین شکل به جای یک مثال منفرد، یک نقشه کمی از حساسیت کل workflow به read depth میدهد.`,
  28: `Supplemental Figure 12 همان نوع خلاصه را برای cell downsampling نشان میدهد. محور افقی fraction سلولهای نگه داشته شده است و هر پنل میپرسد چه زمانی dataset کم شده به اندازه تفاوت Seurat و Scanpy از dataset کامل فاصله میگیرد. metricها مشابه شکل قبلی هستند، از filtering و HVG تا PCA، graph، cluster، UMAP و DE.

این شکل نشان میدهد کم کردن cellها معمولا نسبت به کم کردن readها اثر متفاوتی دارد و برای حفظ agreement باید fraction بیشتری از سلولها نگه داشته شود. پیام عملی این است که sample size در single-cell فقط یک عدد جانبی نیست؛ کاهش تعداد سلول میتواند همان خروجیهایی را جابه جا کند که بعدا برای تفسیر زیستی استفاده میشوند.`,
  29: `Supplemental Figure 13 Cell Ranger را وارد داستان reproducibility میکند. مقایسه بین نسخه 7 و نسخه 6 انجام شده و سپس اثر این تفاوت upstream در کل workflow دنبال میشود: filtering، HVG selection، PCA، graph construction، clustering، UMAP و marker genes. نکته مهم مقاله این است که Cell Ranger v7 به طور پیش فرض intron counts را وارد میکند، در حالی که v6 این کار را نمیکرد.

برای ارائه، این شکل باید یادآوری کند که reproducibility از Seurat یا Scanpy شروع نمیشود. قبل از آن، count matrix توسط ابزار دیگری ساخته شده است و اگر آن ابزار یا defaultهایش عوض شود، ماتریس ورودی عوض میشود. بنابراین حتی اگر downstream analysis ثابت باشد، upstream matrix generation میتواند نتیجه نهایی را تغییر دهد.`,
  30: `Supplemental Figure 14 یک نگاه QC کوتاه به ورودیها میدهد. پنل A توزیع UMI و feature را بین نسخه های Cell Ranger مقایسه میکند و پنل B همین نوع توزیعها را بین داده کامل و read-downsampled نشان میدهد. این plotها قبل از filtering و feature selection قرار دارند و نشان میدهند که landscape ورودی از همان ابتدا میتواند جابه جا شود.

نکته مهم این است که تغییر در UMI یا feature distribution فقط یک تفاوت توصیفی نیست. filtering thresholdها، انتخاب ژنهای متغیر و حتی ساختار downstream ممکن است به این توزیعهای اولیه حساس باشند. بنابراین این شکل پشت صحنه Figure 13 و downsampling را توضیح میدهد: وقتی ورودی تغییر کند، همه مراحل بعدی ممکن است تغییر را حمل کنند.`,
  31: `Supplemental Figure 15 تحلیل downstream زیستی را برای clusterهای 0 تا 9 گسترش میدهد. هر ردیف یک cluster است و ستونها به ترتیب volcano plot Seurat، volcano plot Scanpy، overlap ژنهای biologically significant و agreement در Enrichr pathwayها را نشان میدهند. رنگ قرمز ژنهای upregulated معنی دار، آبی downregulated معنی دار، و خاکستری غیر معنی دار است.

اهمیت این اسلاید این است که مثال cluster 15 در Figure 4 را از حالت یک نمونه انتخابی خارج میکند. اینجا چندین cluster دیده میشوند و در بسیاری از آنها تفاوت بین Seurat و Scanpy در gene calls و pathwayها تکرار میشود. بنابراین پیام زیستی مقاله فقط به یک cluster خاص وابسته نیست؛ defaultهای workflow میتوانند در سطح وسیع تری تفسیر clusterها را تغییر دهند.`,
  32: `Supplemental Figure 16 ادامه همان تحلیل برای clusterهای 10 تا 18 است. چون layout دقیقا مشابه شکل قبلی است، میتوان سریع scan کرد که کجا volcano plotها، overlap ژنهای معنی دار، یا pathway agreement بین Seurat و Scanpy فرق میکند. هدف این نیست که تک تک ژنهای هر ردیف حفظ شوند، بلکه باید الگوی کلی دیده شود.

این اسلاید coverage تحلیل را کامل میکند. همراه Supplementary Figure 15 نشان میدهد که اختلاف downstream بین دو package در چندین cluster دیده میشود. بنابراین وقتی در Figure 4 یک cluster خاص با جزئیات نشان داده شد، این supplement ثابت میکند آن مثال فقط برای نمایش انتخاب نشده، بلکه بخشی از یک الگوی گسترده تر است.`,
  33: `Supplemental Figure 17 اثر random seed و stochasticity را بررسی میکند. نیمه بالا اجزای تصادفی در Seurat را نشان میدهد، مثل neighbor search، Louvain clustering، uwot UMAP و گرافهای UMAP-derived. نیمه پایین اجزای مشابه در Scanpy را بررسی میکند، از جمله NNDescent، Leiden و UMAP-learn. هر metric با تفاوت package-level مقایسه میشود تا اندازه اثر seed قابل تفسیر باشد.

این شکل برای بحث منصفانه خیلی مهم است. ممکن است کسی بگوید اختلاف Seurat و Scanpy فقط از random seed آمده است. این supplement نشان میدهد seed واقعا میتواند خروجی را جابه جا کند، اما همه اختلافهای package-level را توضیح نمیدهد. پس stochasticity یک منبع variability است، ولی علت کامل تفاوتهای اصلی مقاله نیست.`,
  34: `Supplemental Figure 18 یک schematic برای فهمیدن اختلاف KNN و SNN graph است. در پنل A نشان داده میشود که approximate KNN ممکن است بعضی edgeها را اضافه یا از دست بدهد. در پنل B ایده pruning در graph سبک Seurat توضیح داده شده است؛ shared-neighbor structure تعیین میکند کدام edgeها باقی بمانند، و این رفتار برای nodeهای hub و peripheral میتواند متفاوت باشد. Scanpy بیشتر به undirected KNN graph نزدیک است.

این شکل کمک میکند مرحله graph construction را از حالت جعبه سیاه خارج کنیم. حتی اگر PCA دو روش شبیه به نظر برسد، تبدیل آن به graph میتواند متفاوت باشد. چون clustering و گاهی UMAP روی همین graphها ساخته میشوند، تفاوت در edgeها و pruning میتواند مستقیم به تفاوت در cluster structure برسد.`,
  35: `Supplemental Figure 19 مسیر دقیق graph construction را در Seurat و Scanpy کنار هم میگذارد. مسیر بالا نشان میدهد Seurat چگونه KNN و SNN graph را میسازد، prune میکند، ذخیره میکند و برای clustering یا UMAP استفاده میکند. مسیر پایین route مربوط به Scanpy و objectهای گرافی آن را نشان میدهد. متنهای bold نام function یا method هستند و متنهای italic packageهای خارجی را مشخص میکنند.

پیام این diagram این است که اختلاف گراف فقط یک parameter ساده نیست. دو package objectها را در جاهای متفاوتی میسازند، با conventionهای متفاوتی نگه میدارند و در downstream functionهای متفاوت استفاده میکنند. برای همین matching کامل behavior گراف بین Seurat و Scanpy دشوار است و باید architecture workflow را فهمید.`,
  36: `Supplemental Figure 20 بررسی میکند آیا نتیجه اصلی فقط مخصوص dataset PBMC بوده یا در یک سیستم زیستی دیگر هم دیده میشود. اینجا همان مقایسه default Seurat در برابر Scanpy در mouse brain تکرار شده است. پنلها همان ترتیب Figure 1 را دارند: filtering و HVG، PCA، graph neighborhoods، clustering و UMAP، و در آخر DE.

برای ارائه، این اسلاید یک robustness check است. اگر اختلاف فقط در PBMC دیده میشد، ممکن بود به ترکیب خاص آن dataset مربوط باشد. اما در mouse brain هم workflow divergence دیده میشود. بنابراین مقاله ادعا میکند مشکل فقط یک artifact dataset-specific نیست، بلکه میتواند در سیستمهای زیستی متفاوت هم ظاهر شود.`,
  37: `Supplemental Figure 21 همان کنترل external dataset را در non-small cell lung cancer انجام میدهد. دوباره ساختار پنلها مثل Figure 1 است تا بتوان همان checkpointها را دنبال کرد: از filtering و HVG تا PCA، گراف، clustering، UMAP و differential expression. با وجود زمینه بیماری و ترکیب سلولی متفاوت، defaultهای Seurat و Scanpy همچنان خروجیهای قابل اندازه گیری متفاوت میسازند.

این اسلاید ادعای مقاله را گسترده تر میکند. تفاوت packageها فقط در یک dataset آموزشی یا PBMC ساده دیده نشده، بلکه در dataset سرطان ریه هم باقی است. بنابراین برای هر تحلیل single-cell واقعی، مخصوصا تحلیلهایی که قرار است تفسیر زیستی یا clinical داشته باشند، مستند کردن package defaults و نسخه ها اهمیت دارد.`,
  38: `Supplemental Figure 22 variability ناشی از bootstrap را خلاصه میکند. هر پنل یک metric consistency را روی bootstrap datasetهای تولید شده با seedهای مختلف نشان میدهد. رنگ نارنجی Seurat و آبی Scanpy است و خط dashed مثل شکلهای قبلی تفاوت کامل Seurat در برابر Scanpy را نشان میدهد. بنابراین bootstrap variation هم با همان benchmark package-choice مقایسه میشود.

این شکل میپرسد resampling noise داخل یک workflow چقدر بزرگ است. اگر bootstrap variation به خط dashed نزدیک شود، یعنی تغییرات ناشی از نمونه گیری هم اندازه تغییر package است. اگر پایین تر بماند، package choice اثر بزرگ تری دارد. این مقایسه کمک میکند اختلافهای مقاله را در کنار نویز طبیعی نمونه گیری تفسیر کنیم.`,
  39: `Supplemental Figure 23 یک نمونه ملموس از bootstrap را نشان میدهد: یک اجرای bootstrap در Seurat در برابر dataset اصلی. پنلها دوباره مسیر کامل pipeline را طی میکنند، از filtering و HVG تا PCA، graph، clustering، UMAP و DE. این شکل باعث میشود خلاصه های خطی bootstrap در شکل قبلی کمتر انتزاعی به نظر برسند.

برای ارائه، باید گفت این اسلاید نشان میدهد resampling واقعا میتواند خروجی workflow را تکان دهد، اما مقاله آن را در برابر benchmark بزرگ تر Seurat-versus-Scanpy میسنجد. بنابراین هدف این نیست که bootstrap را بی اهمیت نشان دهد؛ هدف این است که اندازه این منبع variability را در کنار اثر package choice قرار دهد.`,
  40: `Table S1 یک جدول عملی برای تنظیمات Seurat است. ردیفها مراحل workflow را نشان میدهند: filtering، normalization، انتخاب HVG، scaling، PCA، SNN، clustering، UMAP و differential expression. ستونها مشخص میکنند default Seurat چیست، با چه argumentهایی میتوان آن را به Scanpy نزدیک کرد، و کجاها behavior دو ابزار compatible، partly incompatible یا incompatible است.

این جدول برای بخش discussion یا methods خیلی مفید است. پیام آن این است که reproducibility فقط با نوشتن «ما از Seurat استفاده کردیم» حاصل نمیشود. باید function، version، default و argumentهای تغییر یافته برای هر مرحله ثبت شوند. این جدول عملا چک لیست Seurat برای کسی است که میخواهد تحلیلش با Scanpy قابل مقایسه یا قابل بازتولید باشد.`,
  41: `Table S2 همان چک لیست را از سمت Scanpy ارائه میکند. ردیفها همان مرحله های workflow هستند و جدول نشان میدهد Scanpy از چه functionها و defaultهایی استفاده میکند، کدام بخشها با Seurat از ابتدا match هستند، کجا با argument change نزدیک میشوند، و کجا اصلا معادل تمیز وجود ندارد. این جدول complement طبیعی Table S1 است.

برای کسی که مقاله را نخوانده، نکته اصلی این است که Scanpy و Seurat فقط دو implementation با اسم متفاوت نیستند. در بعضی مرحله ها مفاهیم، objectها یا defaultها کاملا هم ارز نیستند. بنابراین یک methods section خوب باید نسخه، function call، default و تغییرات argument را دقیق بنویسد، مخصوصا در مرحله هایی مثل graph construction، UMAP و differential expression.`,
};

function noteText(slide) {
  return String(slide.speakerNotes?.text || slide.speakerNotes?.getText?.() || "");
}

function stripPersianSection(note) {
  return String(note || "").replace(/\n\s*Persian explanation:\s*[\s\S]*$/i, "").trim();
}

function hasRequiredEnglishSections(note) {
  return /Ready-to-read note:/i.test(note) && /\n\s*Legend:/i.test(note) && /\n\s*Source:/i.test(note);
}

function persianWordCount(text) {
  return String(text || "").split(/\s+/).filter(Boolean).length;
}

function audit(presentation) {
  const missing = [];
  const tooShort = [];
  const noPersianLetters = [];
  const missingEnglish = [];
  for (const slideNumber of Object.keys(FARSI).map(Number).sort((a, b) => a - b)) {
    const slide = presentation.slides.items[slideNumber - 1];
    const note = noteText(slide);
    const match = note.match(/Persian explanation:\s*([\s\S]*)$/i);
    if (!hasRequiredEnglishSections(note)) {
      missingEnglish.push(slideNumber);
    }
    if (!match || !match[1].trim()) {
      missing.push(slideNumber);
      continue;
    }
    const text = match[1].trim();
    if (!/[\u0600-\u06FF]/.test(text)) {
      noPersianLetters.push(slideNumber);
    }
    if (persianWordCount(text) < 80) {
      tooShort.push({ slide: slideNumber, words: persianWordCount(text) });
    }
  }
  return { missing, tooShort, noPersianLetters, missingEnglish };
}

const presentation = await PresentationFile.importPptx(await FileBlob.load(SOURCE));

for (const [slideNumber, explanation] of Object.entries(FARSI)) {
  const slide = presentation.slides.items[Number(slideNumber) - 1];
  if (!slide) {
    throw new Error(`Missing slide ${slideNumber}`);
  }
  const base = stripPersianSection(noteText(slide));
  if (!hasRequiredEnglishSections(base)) {
    throw new Error(`Slide ${slideNumber} does not have expected English note sections.`);
  }
  slide.speakerNotes.setText(`${base}\n\nPersian explanation:\n${explanation.trim()}`);
}

const preExportAudit = audit(presentation);
if (
  preExportAudit.missing.length ||
  preExportAudit.tooShort.length ||
  preExportAudit.noPersianLetters.length ||
  preExportAudit.missingEnglish.length
) {
  throw new Error(`Pre-export Persian note audit failed: ${JSON.stringify(preExportAudit, null, 2)}`);
}

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(OUT);

const reloaded = await PresentationFile.importPptx(await FileBlob.load(OUT));
const postExportAudit = audit(reloaded);
if (
  postExportAudit.missing.length ||
  postExportAudit.tooShort.length ||
  postExportAudit.noPersianLetters.length ||
  postExportAudit.missingEnglish.length
) {
  throw new Error(`Post-export Persian note audit failed: ${JSON.stringify(postExportAudit, null, 2)}`);
}

console.log(
  JSON.stringify(
    {
      source: SOURCE,
      output: OUT,
      slideCount: reloaded.slides.count,
      updatedSlides: Object.keys(FARSI).map(Number).sort((a, b) => a - b),
      postExportAudit,
    },
    null,
    2,
  ),
);

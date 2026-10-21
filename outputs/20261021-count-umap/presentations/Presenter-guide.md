# From Count Matrix to UMAP: What Happens in Between?

Main route: slides 1–22, 25 minutes. Discussion: slide 23, 5 minutes. Slides 24 onward are backup material.

Scientific anchor: Luecken & Theis (2019), doi:10.15252/msb.20188746.

Visible slide text is English. Every slide has an English script and Persian preparation notes.

## 1. From Count Matrix to UMAP: What Happens in Between?

30 seconds

A UMAP plot is the end of a chain of decisions. Today we will follow that chain from molecular counts to a map of cells, using the tutorial by Luecken and Theis as our main paper. The central question is why we usually reduce dimensions twice. PCA creates a compact representation of expression variation. UMAP uses neighborhood relationships in that representation to place cells on a page. We will spend most of our time on the connection between those steps, then discuss which conclusions the final picture can support.

هدف این جلسه فهم منطق مسیر است: دادهٔ شمارش چگونه به شباهت بین سلول‌ها و در نهایت به یک تصویر تبدیل می‌شود؟ مقالهٔ Luecken و Theis پایهٔ بحث است. تفاوت بین خلاصه‌سازی برای تحلیل و کاهش بُعد برای نمایش، نخ اصلی ارائه است.

## 2. The paper separates preprocessing from biological inference

60 seconds

Figure 1 places our question inside the full analysis workflow. On the left, counts pass through quality control, normalization, correction when needed, feature selection and visualization. On the right are biological tasks such as annotation and trajectory inference. Today we stop at the representation that makes these tasks possible. The examples come from mouse intestinal epithelium, so the same dataset will connect the stages. This is a review and tutorial, rather than a benchmark proving one universal pipeline. Its most useful distinction for today is between preparing a representation and drawing biological conclusions from it.

شکل ۱ نقشهٔ کل مقاله است. امروز بخش چپ شکل و منطق رسیدن به UMAP را باز می‌کنیم. سمت راست یادآوری می‌کند که دیدن خوشه تازه شروع تفسیر زیستی است. این مقاله مرور آموزشی است و برندهٔ قطعی همهٔ روش‌ها را تعیین نمی‌کند.

## 3. A cell begins as a vector of molecular counts

60 seconds

Each row here is one cell and each gene is one coordinate. The table is intentionally tiny, but the actual matrix may contain thousands of cells and many thousands of genes. A count measures detected molecules in a UMI protocol, rather than a direct concentration. A zero can reflect lack of expression or limited sampling. Cells A and B have the same relative pattern in the displayed genes, but B has twice the library size. If we compare raw counts directly, sequencing depth can therefore influence which cells appear similar before we have learned anything about their identity.

هر سلول یک بردار در فضای ژن‌هاست. مثال عمداً ساختگی است. A و B از نظر نسبت بیان ژن‌های نشان‌داده‌شده مشابه‌اند ولی عمق شمارش B دو برابر است. صفر الزاماً اثبات خاموش‌بودن ژن نیست. جهت ماتریس در نرم‌افزارها ممکن است برعکس باشد، اما مفهوم یکسان است.

## 4. QC changes which cells enter the geometry

75 seconds

Figure 2 shows the first decision that changes our eventual map: which barcodes we retain. Panels A and B show count and gene distributions, and panel C shows the ranked count curve. Panel D combines detected genes, count depth and mitochondrial fraction. The red cutoffs belong to this particular intestinal dataset. They are examples, rather than universal thresholds. A high count can reflect a large cell or a doublet, and a high mitochondrial fraction needs tissue context. Removing a cell here removes its contribution to every later covariance, principal component and neighborhood relationship.

در QC فقط داده را تمیز نمی‌کنیم؛ جمعیت سلول‌هایی را تعیین می‌کنیم که هندسهٔ بعدی بر اساس آن‌ها ساخته می‌شود. آستانه‌های قرمز شکل مخصوص همین داده‌اند. کم‌بودن RNA می‌تواند ویژگی یک جمعیت واقعی باشد، بنابراین تعداد ژن، UMI و درصد میتوکندری را با هم و در زمینهٔ هر نمونه می‌خوانیم.

## 5. Normalization changes the meaning of a distance

75 seconds

Library-size normalization divides a gene count by the total count in that cell, then applies a common scale factor. A and B now both have a normalized value of forty for this gene, and log one plus forty is approximately three point seven one four. The log transformation compresses the contribution of very large values. This changes the geometry on which PCA operates. The procedure assumes that relative expression is a useful comparison, and it does not recover absolute RNA content. Model-based alternatives, including SCTransform or residual-based workflows, change the normalization route but still require us to understand what representation enters dimensionality reduction.

این مثال عددی نشان می‌دهد نرمال‌سازی چه چیزی را برابر می‌کند: سهم نسبی یک ژن در کتابخانه. لگاریتم دامنهٔ اعداد بزرگ را فشرده می‌کند. بیان نسبی با تعداد مطلق RNA برابر نیست. مسیر کلاسیک log-normalization را برای آموزش دنبال می‌کنیم؛ SCTransform یک مسیر جایگزین است.

## 6. HVG selection chooses genes, scaling changes their weight

75 seconds

Feature selection and scaling do different jobs. HVG selection chooses which gene columns will define our representation, ideally emphasizing variability beyond the expected mean–variance relationship. Scaling changes the relative contribution of the selected genes. Centering removes a gene-specific mean, while unit-variance scaling expresses deviations in standard-deviation units. Scaling can prevent highly expressed genes from dominating, but it can also amplify noisy low-expression features. The original count matrix should remain available. The order and input depend on the HVG method, because some methods expect counts and others expect log-normalized expression. The choices here determine what variance PCA will regard as important.

انتخاب HVG یعنی انتخاب محورهای ژنیِ ورودی؛ scaling یعنی تغییر وزن همان محورها. حذف ژن از محاسبهٔ هندسه نباید به معنی از دست‌دادن شمارش اصلی برای تحلیل بعدی باشد. روش‌های HVG ورودی یکسانی ندارند و ترتیب نسخهٔ کلاسیک را نباید به همهٔ روش‌ها تعمیم داد.

## 7. Batch correction can change the neighborhood structure

60 seconds

Figure 3 compares a UMAP before and after batch correction, with cells colored by sample. The sample-specific arrangement changes after correction. This demonstrates that our map depends on how we define the representation before visualization. It does not prove that the corrected biology is right. We still need to check whether comparable cell identities align and whether real condition-specific states remain distinct. When condition and batch are completely confounded, a computational correction cannot establish which difference is technical. In current workflows, integration may adjust expression, a latent representation or a graph, so it does not have one mandatory position in every pipeline.

رنگ‌ها نمونه‌ها را نشان می‌دهند. مخلوط‌شدن رنگ‌ها بعد از correction به‌تنهایی معیار موفقیت نیست؛ باید هویت‌ها و حالت‌های زیستی حفظ شوند. اگر گروه زیستی با batch کاملاً یکی باشد، الگوریتم نمی‌تواند منشأ اختلاف را اثبات کند. جای integration به نوع روش وابسته است.

## 8. Many genes share a smaller number of variation patterns

60 seconds

A cell has thousands of gene measurements, but many of those measurements vary together. A coordinated expression program can therefore occupy fewer effective dimensions than its number of genes. We also observe technical variation and measurement noise. PCA combines correlated features into a compact set of coordinates and can reduce the influence of some noise when we truncate the representation. The numbers on this slide describe an illustrative route, not fixed settings. Retaining thirty PCs is still a thirty-dimensional analysis. The final move to two UMAP coordinates serves the separate goal of making the relationships visible on a screen.

تعداد ژن با تعداد الگوهای مستقل تغییرات یکی نیست. یک برنامهٔ زیستی ممکن است ده‌ها ژن را هم‌جهت تغییر دهد. PCA از هم‌بستگی‌ها برای خلاصه‌کردن استفاده می‌کند. عدد ۳۰ فقط مثال است؛ قبل از UMAP هنوز در یک فضای ۳۰بعدی مقایسه می‌کنیم، نه روی دو محور اول PCA.

## 9. PCA rotates the axes, truncation removes directions

90 seconds

In the left plot, the two genes vary together, so most of the spread lies along a diagonal direction. PCA defines PC1 along the direction with the greatest variance and PC2 along the perpendicular direction with the next greatest variance. The right plot expresses the cloud in these new coordinates. Keeping both components is a rotation and preserves the full information in this two-dimensional example. Dimensionality reduction occurs when we retain PC1 and discard the smaller perpendicular variation. The orange ticks illustrate that projection. In real data, discarded variation can include noise and useful biology, which is why truncation needs inspection.

خود چرخشِ همهٔ محورها چیزی را حذف نمی‌کند. کاهش بُعد وقتی اتفاق می‌افتد که فقط چند PC را نگه داریم. PC1 بیشترین واریانس را می‌گیرد و PC2 در جهت عمود، بیشترین واریانس باقیمانده را. مقدار کمِ واریانس معادل بی‌ارزش‌بودن زیستی نیست؛ یک جمعیت نادر ممکن است در مؤلفه‌های بعدی دیده شود.

## 10. Loadings describe genes, scores locate cells

90 seconds

The matrix equation separates two objects that are often confused. Z contains the centered, and possibly scaled, expression values. The columns of V are loading vectors, with one weight per gene. Multiplying them gives T, the scores that locate each cell in PC space. In the simple example, a cell high in both genes receives a high PC1 score. A loading can also be negative, so a component may contrast two expression programs. A principal component is therefore a weighted combination of genes rather than one selected gene or an assigned cell type. Loading inspection connects the numerical direction back to possible biology.

loading وزن ژن در یک مؤلفه است و score مختصات یک سلول روی همان مؤلفه. از روی ژن‌های با وزن مثبت و منفی می‌توان معنای احتمالی PC را بررسی کرد. PC معادل یک ژن یا یک خوشه نیست. تغییر علامت کل یک PC هیچ معنای زیستی تازه‌ای ایجاد نمی‌کند.

## 11. Truncating PCA changes the distances between cells

75 seconds

The retained scores define the distances used for neighbor search. With an orthonormal PCA basis and every available component retained, PCA preserves Euclidean distances in the centered input space. Once we truncate, only differences along the retained components contribute. This explains both the benefit and the risk. If discarded directions mainly contain noise, neighborhoods can become more useful. If a rare cell state differs along a discarded direction, its separation can disappear. Standard PCA scores retain their variance differences. Whitening would reweight the components and change the geometry again, so it should not be confused with ordinary PCA.

فاصلهٔ دو سلول در فضای PCA از اختلاف scoreها در مؤلفه‌های نگه‌داشته‌شده محاسبه می‌شود. با نگه‌داشتن همهٔ مؤلفه‌ها، چرخش فاصله را عوض نمی‌کند؛ با حذف مؤلفه‌ها، سهم آن جهت‌ها حذف می‌شود. بنابراین PCA تضمین پاک‌کردن نویز نیست و ممکن است بخشی از سیگنال را هم حذف کند.

## 12. The elbow is a starting point for choosing PCs

75 seconds

Panel F of Figure 4 ranks the variance explained by the first thirty-one principal components. The steep decline becomes flatter around PCs five to seven in this intestinal dataset. The paper uses that elbow as a practical heuristic. It does not establish a universal PC count. A better decision also considers the genes loading on the components, associations with depth or batch, and the stability of relevant neighborhoods over plausible choices. A low-variance component may describe a small population. Conversely, a dominant component may reflect technical variation. The goal is a useful representation for the biological question, rather than maximizing retained variance alone.

آرنج نمودار برای این داده حدود PCهای ۵ تا ۷ است. این عدد را برای دادهٔ خودمان کپی نمی‌کنیم. loadings، رابطه با QC و batch، و پایداری همسایگی‌ها را در چند انتخاب معقول بررسی می‌کنیم. بیشترین واریانس الزاماً مهم‌ترین زیست‌شناسی نیست.

## 13. A strong PC can reflect sequencing depth

60 seconds

Panel A shows the first two PCs of the same intestinal cells. The original figure colors cells by count depth, and high-depth cells occupy a distinct part of the representation. PCA has no knowledge of cell types or technical quality. It simply finds variation in the matrix we supplied. A depth association can indicate residual technical structure, but depth can also correlate with real cell size or activity. This is a reason to inspect the relationship and the relevant genes, rather than automatically delete PC1. The question is whether the component captures the biology we want to compare.

در شکل اصلی، رنگ نقاط عمق شمارش است. PCA برچسب زیستی نمی‌شناسد و حتی PC1 می‌تواند با عمق یا batch همراه باشد. این مشاهده دلیل بررسی است، نه مجوز حذف خودکار PC1. عمق گاهی با اندازه یا فعالیت واقعی سلول مرتبط است.

## 14. Neighbor search turns PC coordinates into local relationships

75 seconds

The next object is a neighborhood graph. In this toy example, A is the reference cell and B and C are its two nearest neighbors in three-dimensional PC space. In a real analysis, the same search uses all retained PCs, not just the first two on a plot. Repeating the search for each cell produces local relationships. UMAP turns local distances into weighted, fuzzy connections, adapting to local scales. This is why the number of neighbors matters. A graph emphasizes who is connected to whom, and does not retain every pairwise distance with equal importance.

گام واسطِ فراموش‌شده بین PCA و UMAP ساخت رابطهٔ همسایگی است. جدول مثال ساختگی سه‌بعدی است و سلول A همسایهٔ خودش حساب نمی‌شود. در دادهٔ واقعی از همهٔ PCهای انتخاب‌شده استفاده می‌شود. وزن یال در UMAP صرفاً خود فاصلهٔ خام نیست و به مقیاس محلی هم بستگی دارد.

## 15. UMAP fits a two-dimensional layout to a weighted graph

90 seconds

UMAP first builds a weighted description of local relationships in the input space. It then searches for low-dimensional coordinates whose connection strengths resemble that description. Strong connections create attraction in the optimization, while negative sampling approximates a repulsive contribution for other pairs. A flat map cannot satisfy all relationships in a complex dataset, so the algorithm makes tradeoffs. Initialization, the optimization and the input graph all affect the arrangement. The axes are coordinates of this fitted layout, and do not have PCA-like gene loadings. UMAP can embed into more dimensions, but here we use two because our purpose is visualization.

UMAP ابتدا رابطه‌های محلی وزن‌دار می‌سازد و سپس مختصاتی می‌یابد که آن رابطه‌ها را در بُعد کم تقریب بزند. نمایش تخت ناچار به مصالحه است. محور UMAP1 وزن ژنی قابل‌خواندن مثل یک PC ندارد. UMAP الزاماً محدود به دو بُعد نیست، ولی اینجا دو بُعد را برای نمایش می‌خواهیم.

## 16. The same cells produce different views of structure

75 seconds

Figure 4 compares several views of the same cells, with count depth supplying the shared color scale. Panel A uses the first two PCs. Panels B through E show t-SNE, diffusion maps, UMAP and a force-directed layout. Panel F gives the PCA variance curve. The shapes differ because the methods optimize different properties. A compact PCA representation can still contain useful structure even when its first two coordinates overlap visually. UMAP can make local relationships easier to see in two dimensions, but the appearance does not independently validate the cell states or guarantee faithful distances between separated groups.

همان سلول‌ها و همان متغیر رنگ در چند روش، شکل‌های متفاوتی می‌سازند. نامرتب‌بودن دو PC اول به معنی بی‌فایده‌بودن فضای چندده‌بعدی PCA نیست. UMAP برای نمایش دوبعدی مناسب‌تر است، اما تصویر جذاب‌تر صحت زیستی یا فاصلهٔ دقیق بین جزیره‌ها را ثابت نمی‌کند.

## 17. min_dist changes how tightly points pack

60 seconds

These two official tutorial plots use the same synthetic color dataset and change the minimum-distance setting. The smaller value permits more tightly packed groups, while the larger value spreads points more broadly. The example is deliberately not a cell dataset, which helps separate an algorithmic effect from a biological interpretation. In a workflow with a fixed input graph, changing min_dist changes the display rather than the expression matrix. The number of neighbors has a different role because it changes the local relationships considered. A more compact island should therefore prompt a parameter check before a claim about a sharper cell identity.

هر دو تصویر از یک دادهٔ مصنوعی رنگ‌ها هستند، نه سلول واقعی. min_dist کوچک‌تر اجازهٔ تجمع فشرده‌تر می‌دهد. با گراف ثابت، این تنظیم ظاهر نقشه را تغییر می‌دهد. n_neighbors نقش دیگری دارد و تعریف همسایگی را تغییر می‌دهد؛ این دو پارامتر قابل‌جایگزینی نیستند.

## 18. PCA and UMAP solve different problems

75 seconds

The sequence makes sense because the objectives differ. PCA provides a compact linear representation with interpretable gene loadings, and UMAP makes a local relationship structure visible in two dimensions. Applying UMAP directly to all genes is possible, so PCA is not a mathematical prerequisite. In typical single-cell workflows, PCA offers computational efficiency and a way to reduce redundant or noisy directions before neighbor search. These benefits depend on the preprocessing and the retained components. UMAP cannot restore a rare-state signal that we removed earlier. Two reductions are useful when their different objectives match our analysis and communication needs.

PCA پیش‌شرط ریاضی UMAP نیست. استفادهٔ رایج از این زنجیره به دلیل تفاوت هدف‌هاست: ساخت یک نمایش فشرده و بعد قابل‌دیدن‌کردن روابط آن. اگر در HVG یا حذف PCها سیگنالی را دور ریخته باشیم، UMAP آن را برنمی‌گرداند.

## 19. Clustering and UMAP branch from related representations

60 seconds

A common misunderstanding is that the visible UMAP islands are the input to routine single-cell clustering. In typical workflows, graph clustering and visualization branch from related high-dimensional representations. Scanpy can use a neighborhood graph for both Leiden clustering and UMAP. Seurat commonly uses an SNN graph for clustering while RunUMAP builds its own neighborhoods from the specified PCs. These graphs need not be identical. Once labels have been assigned, we color the UMAP by those labels. Changing the layout alone can therefore alter the picture without changing the labels. Other clustering strategies exist, but this is the standard distinction relevant here.

معمولاً خوشه‌بندی روی جزیره‌های دوبعدی UMAP انجام نمی‌شود. در Scanpy گراف همسایگی می‌تواند برای هر دو استفاده شود. در Seurat گراف SNN خوشه‌بندی و همسایگی داخلی RunUMAP الزاماً یکی نیستند. رنگ‌کردن نقاط با برچسب خوشه یعنی نتیجهٔ تحلیل دیگری را روی نقشه نشان می‌دهیم.

## 20. UMAP distances need evidence beyond the picture

60 seconds

The map is useful for generating questions. Nearby cells may share local expression structure, but we should check the relevant markers and representation. A gap can depend on preprocessing, neighborhood choices and the layout. The area of an island does not estimate the fraction of cells in that population, because UMAP does not preserve density in a straightforward way. An apparent bridge does not establish direction, developmental time or lineage. Strong biological claims need independent evidence and a design that supports them. A useful journal-club question is whether the conclusion survives plausible choices of features, PCs and neighborhood settings.

از شکل جزیره نباید فراوانی جمعیت، مرز قطعی نوع سلولی یا مسیر زمانی استنتاج کنیم. اندازهٔ جزیره با تعداد نقاط یکی نیست و UMAP چگالی را ساده و مستقیم حفظ نمی‌کند. نتیجهٔ زیستی را با مارکرها، طراحی مطالعه، تکرارهای زیستی و پایداری نمایش بررسی می‌کنیم.

## 21. Different analysis tasks use different data layers

45 seconds

Table 1 in the paper reminds us that analysis produces multiple data layers rather than a single progressively better matrix. We preserve raw counts because library-aware models and sample-level pseudobulk analyses need them. We use appropriate gene-expression values to display markers, and PC scores or another latent representation to compare cells. UMAP coordinates are mainly a visualization product in this workflow. Differential-expression methods have different input requirements, and a biological comparison also requires the correct experimental unit. The practical lesson is to retain provenance and choose the representation for the task instead of feeding the latest output into every subsequent analysis.

این جدول برداشت آموزشی از Table 1 است، نه بازنویسی لفظی همهٔ توصیه‌های سال ۲۰۱۹. یک ماتریسِ نهایی برای همهٔ کارها نداریم. شمارش خام، بیان نرمال‌شده، PCA و UMAP هر کدام کاربرد خود را دارند. برای مقایسهٔ شرایط، واحد آزمایشی و تکرار زیستی هم ضروری‌اند.

## 22. The full chain explains the final map

75 seconds

We can now answer the title question. Counts record sampled molecules. Quality control chooses the cells we analyze, and normalization defines a useful comparison scale. Feature selection and scaling decide which gene differences matter. PCA summarizes that representation, and neighbor search expresses local similarity. UMAP fits display coordinates to those relationships. This order is a useful classical route, with method-specific alternatives at several steps. The final map inherits each upstream decision. Understanding the sequence allows us to diagnose a surprising plot by returning to the representation that produced it, rather than adjusting the appearance alone.

زنجیره را یک بار با زبان ساده جمع می‌کنیم. هر مرحله نوعی تصمیم دربارهٔ مشاهده، مقایسه یا نمایش است. این ترتیب مسیر کلاسیک آموزشی است و بسته به روش می‌تواند تغییر کند. برای عیب‌یابی یک UMAP عجیب باید به ورودی و همسایگی برگردیم، نه فقط تنظیم ظاهر.

## 23. Discussion: which decisions would change your conclusion?

300 seconds

The paper gives us a strong organizing framework, especially the distinction between summarization and visualization. Its limitations are equally useful for discussion. It is a 2019 review, its worked example is one tissue system, and a heuristic such as an elbow does not guarantee biological completeness. Consider these three observations. A rare population vanishes when we retain fewer PCs. A cluster becomes tighter after changing min_dist. Integration mixes samples but removes a disease-associated state. In each case, which object changed, what evidence would we inspect next, and what conclusion can we support without treating the UMAP picture as ground truth?

برای پنج دقیقهٔ آخر، سه سناریو را مقایسه می‌کنیم. ناپدیدشدن جمعیت با کاهش PC مسئلهٔ نمایشِ تحلیلی است؛ فشرده‌شدن با min_dist بیشتر مسئلهٔ چیدمان است؛ حذف حالت بیماری بعد از integration نیازمند ارزیابی حفظ زیست‌شناسی است. مقاله چارچوب فکری مفیدی دارد، ولی مرور ۲۰۱۹ و یک مثال بافتی تضمین نسخهٔ عمومی نیستند.

## 24. Appendix: complete figures and implementation

Appendix

The remaining slides retain the downstream figures, expanded-view evidence, the appendix example and implementation references for discussion.

از اینجا اسلایدهای پشتیبان شروع می‌شوند و در زمان اصلی ارائه قرار ندارند. تمام شکل‌های باقی‌ماندهٔ مقاله و ضمائم و دو مثال کد در این بخش‌اند.

## 25. Downstream tasks ask different biological questions

Appendix

Figure 5 groups downstream tasks by the level at which they ask a question. The cell-level branch includes clustering and trajectory analysis, while the gene-level branch includes differential expression, gene sets and regulatory networks. Cluster analysis also includes annotation and comparisons of composition. The same UMAP-like arrangement appears under several tasks, but the picture is only a common display. The appropriate input, statistical unit and interpretation differ across tasks. This figure explains why reaching a visually organized map is an intermediate result. We still need a defined biological question and a method whose assumptions match that question.

شکل ۵ کاربردهای پایین‌دستی را در سطح سلول و ژن جدا می‌کند. استفاده از یک نقشهٔ مشابه برای چند سؤال به معنی یکسان‌بودن ورودی یا روش آماری آن‌ها نیست. خوشه، trajectory، ترکیب جمعیتی و بیان تفاضلی هر کدام سؤال متفاوتی دارند.

## 26. Markers support the interpretation of cluster labels

Appendix

Figure 6 connects graph-derived clusters to biological interpretation. Panel A shows annotated cell identities in a UMAP representation. Panel B overlays example marker genes, including Slc12a2 for stem cells, Arg2 for enterocytes, Tff3 for goblet cells and Defa24 for Paneth cells. Panel C compares the distribution of proximal and distal enterocytes in the map. The clustering supplies groups, while marker expression and prior biological knowledge help interpret those groups. A separated island alone cannot establish a cell type. Expression overlays should come from an appropriate gene-expression layer, and the labels need supporting evidence beyond the visual separation.

شکل ۶ از رنگ خوشه به معنای زیستی آن می‌رسد. مارکرها برای تفسیر برچسب‌ها کمک می‌کنند، نه اینکه جزیرهٔ UMAP به‌تنهایی نوع سلولی را اثبات کند. پنل C محل enterocyteهای proximal و distal را مقایسه می‌کند. پنل‌ها باید همراه با legend و شناخت بافت خوانده شوند.

## 27. Trajectory inference adds assumptions beyond a UMAP layout

Appendix

Figure 7 shows several representations used to study a continuous process. Panel A maps enterocyte relationships, and panel B fits Slingshot trajectories through clusters in PCA space. Panel C relates cluster distributions to pseudotime. Panel D shows a graph abstraction projected onto a UMAP, and panel E displays gene-expression patterns over the inferred ordering. These panels combine a representation with additional modeling assumptions. The trajectory does not emerge merely because a UMAP contains an apparent bridge. Root selection, connectivity and the biological interpretation of pseudotime all matter. Cross-sectional transcriptional similarity alone does not directly measure elapsed time or prove lineage.

شکل ۷ نشان می‌دهد مسیر و pseudotime از مدل و فرض‌های اضافه می‌آیند. Slingshot در پنل B در فضای PCA کار می‌کند و پنل D گراف انتزاع‌شده را روی UMAP نمایش می‌دهد. شباهت رونویسی یا یک پل ظاهری، به‌تنهایی زمان واقعی یا دودمان را اثبات نمی‌کند.

## 28. Published workflows used a range of HVG counts

Appendix

Figure EV1 summarizes a manual survey of published workflows. The horizontal axis gives the logarithm of cell number, the vertical axis gives the number of highly variable genes, and colors indicate technology. The survey shows considerable variation in choices across studies. It provides context for the paper’s suggested range, rather than a controlled test of the optimal number of genes. Dataset size alone does not identify a universally correct setting. Tissue complexity, the signal of interest and the feature-selection method also matter. The underlying survey records are supplied in Dataset EV1.

این شکل پیمایش دستی مقاله‌هاست و benchmark کنترل‌شده نیست. تنوع انتخاب تعداد HVG را نشان می‌دهد. از روی تعداد سلول‌ها به‌تنهایی نمی‌توان بهترین تعداد ژن را تعیین کرد.

## 29. Correction and denoising change gene-level variation

Appendix

Figure EV2 shows changes in the coefficient of variation after correction. Panels A and B use ComBat batch correction for intestinal epithelium and embryonic stem-cell data. Panels C and D use DCA denoising for the corresponding datasets. Negative changes indicate a decrease in the coefficient of variation. The distributions make clear that processing can alter gene-level variation substantially. They do not by themselves establish which changes remove technical noise and which affect useful biology. This matters for our central topic because variance-changing operations upstream can alter the principal components and the neighborhoods used for visualization.

ردیف بالا تغییر CoV بعد از ComBat و ردیف پایین بعد از DCA را نشان می‌دهد. مقادیر منفی یعنی کاهش CoV. تغییر واریانس به‌تنهایی به معنی بهبود نیست؛ باید معلوم شود نویز فنی کم شده یا سیگنال زیستی آسیب دیده است.

## 30. Clustering random data can create apparent marker evidence

Appendix

The appendix provides a useful warning about circular validation. The authors simulated a single population without differentially expressed genes, yet graph clustering produced two groups. Panel A shows those groups in PCA space. Panel B shows an excess of small p-values when marker testing uses the same data that helped define the groups. The test ignores the selection created by clustering. A visually distinct arrangement and significant markers from the same data therefore cannot serve as independent proof of a biological population. External validation or a design that accounts for the selection is needed for stronger claims.

دادهٔ شبیه‌سازی‌شده یک جمعیت دارد و DE واقعی ندارد، اما خوشه‌بندی دو گروه پیدا می‌کند. آزمون مارکر روی همان داده می‌تواند p-value کوچک بدهد. این اعتبارسنجی مستقل نیست و انتخاب ناشی از خوشه‌بندی را نادیده می‌گیرد. نکته برای نقد ادعاهای زیستی از روی نقشه مهم است.

## 31. The paper assigns different uses to each processing stage

Appendix

This editable table retains the stages and application groups listed in Table 1 of the original paper. Measured, corrected and reduced layers serve different purposes. The authors distinguish technical correction from removal of biological covariates, and distinguish feature selection from dimensional summarization. These categories are useful, but should be interpreted in the context of the specific method. A current count model may require raw counts, and a modern integration model may produce a latent representation rather than corrected expression. The enduring lesson is to select the data layer deliberately and retain the information required by each downstream analysis.

در این جدول پوشش کامل مراحل Table 1 مقاله حفظ شده و متن کوتاه شده است. تقسیم‌بندی measured، corrected و reduced مفید است، ولی کاربرد دقیق به فرض روش وابسته است. پیشنهادهای ۲۰۱۹ را به همهٔ مدل‌های جدید تعمیم نمی‌دهیم.

## 32. A classic Seurat route makes the representations explicit

Appendix

This code spells out the classical log-normalized route in Seurat. Variable features define the genes entering PCA, and the explicit dimensions identify the PC coordinates used downstream. FindNeighbors constructs the graph used for clustering. RunUMAP, when given a reduction and dimensions this way, computes its own neighborhood-based embedding from those coordinates. The final labels can then be displayed on the map. The numbers here are illustrative, and the same chosen dimensions should be justified for the dataset. This snippet assumes quality control is complete and preserves the original counts in the object. SCTransform provides an alternative preprocessing route.

کد برای توضیح اتصال مراحل است و روی دادهٔ کاربر اجرا نشده است. QC از قبل فرض شده. ۲۰۰۰ ژن، ۳۰ PC و resolution=0.5 نسخهٔ قطعی نیستند. RunUMAP در این شکل، مستقیماً گراف SNN خوشه‌بندی را دریافت نمی‌کند و از PCA همسایگی می‌سازد.

## 33. A classic Scanpy route exposes the shared neighbor graph

Appendix

This Scanpy example keeps a count layer before normalization, then selects variable genes using a flavor that expects log-normalized expression. PCA uses the explicit highly-variable mask. The neighbor step creates the graph that UMAP and Leiden can both use, which makes the branch visible in the code. Running the visualization before or after assigning cluster labels does not make those labels depend on the two-dimensional layout. Scaling can increase memory use, so this compact teaching example should be adapted for large sparse data. These are illustrative settings rather than a fitted analysis, and the dependencies and input assumptions still require checking.

برای flavor="seurat" ورودی log-normalized است؛ flavor="seurat_v3" انتظار شمارش دارد و جایش فرق می‌کند. گراف neighbors می‌تواند ورودی مشترک UMAP و Leiden باشد. نگه‌داشتن counts برای تحلیل بعدی مهم است. scaling ممکن است حافظهٔ زیادی مصرف کند و این کد مثال آموزشی است.

## 34. Sources and the scope of the update

Appendix

The 2019 tutorial remains the scientific anchor for this presentation. All seven main figures, both expanded-view figures and the appendix figure are retained, with their source legends in the notes. The current documentation clarifies implementation details such as the neighborhood graph and the inputs expected by specific functions. The arithmetic tables and PCA geometry are illustrative teaching additions rather than new experimental results. Recommendations that depend on software or data type should be checked against the documentation and the biological question, while the distinction between summarization, neighborhood construction and visualization remains the central organizing idea.

مرز منابع روشن است: شکل‌ها و چارچوب مفهومی از مقاله‌اند، توضیح کد از مستندات رسمی فعلی و مثال‌های عددی و هندسی ساختهٔ آموزشی‌اند. این ارائه benchmark تازه‌ای انجام نمی‌دهد و نتیجهٔ آزمایش جدیدی ادعا نمی‌کند.

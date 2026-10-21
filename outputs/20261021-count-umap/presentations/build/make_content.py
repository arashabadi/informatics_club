import json,re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
A=ROOT/'assets'
paper=(A/'paper.txt').read_text()
ev=(A/'expanded-view.txt').read_text()
app=(A/'appendix.txt').read_text()
PAPER='Luecken MD & Theis FJ. Molecular Systems Biology 15:e8746 (2019). doi:10.15252/msb.20188746. Figures reproduced under CC BY 4.0.'
SCAN='https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html'
SEURAT='https://satijalab.org/seurat/articles/pbmc3k_tutorial'
UMAP='https://umap-learn.readthedocs.io/en/latest/parameters.html'
ALGO='https://umap-learn.readthedocs.io/en/latest/how_umap_works.html'

def legend(text,label):
    m=re.search(r'(?m)^\s*'+re.escape(label)+r'\.',text)
    if not m:return label
    s=text[m.start():].strip()
    return re.sub(r'\s+',' ',re.split(r'\n\s*\n\s*\n',s,maxsplit=1)[0])

S=[]
def add(title,kind='text',sec=60,**kw):
    s=dict(title=title,kind=kind,seconds=sec,source=PAPER);s.update(kw);S.append(s);return s
def fig(title,num,sec=60,**kw):
    label='Figure '+str(num)
    return add(title,'figure',sec,image=str(A/f'Fig{num}.png'),label=label,legend=legend(paper,label),**kw)

add('From Count Matrix to UMAP:\nWhat Happens in Between?','cover',30,
    subtitle='Why PCA comes before UMAP',
    footer='Arash Bagherabadi\nscTalk 8  /  October 21, 2026',
    talk='A UMAP plot is the end of a chain of decisions. Today we will follow that chain from molecular counts to a map of cells, using the tutorial by Luecken and Theis as our main paper. The central question is why we usually reduce dimensions twice. PCA creates a compact representation of expression variation. UMAP uses neighborhood relationships in that representation to place cells on a page. We will spend most of our time on the connection between those steps, then discuss which conclusions the final picture can support.',
    fa='هدف این جلسه فهم منطق مسیر است: دادهٔ شمارش چگونه به شباهت بین سلول‌ها و در نهایت به یک تصویر تبدیل می‌شود؟ مقالهٔ Luecken و Theis پایهٔ بحث است. تفاوت بین خلاصه‌سازی برای تحلیل و کاهش بُعد برای نمایش، نخ اصلی ارائه است.')

fig('The paper separates preprocessing from biological inference',1,60,
    points=['Our route ends at visualization','The dataset is mouse intestinal epithelium','25 minutes of explanation, 5 minutes of discussion'],
    talk='Figure 1 places our question inside the full analysis workflow. On the left, counts pass through quality control, normalization, correction when needed, feature selection and visualization. On the right are biological tasks such as annotation and trajectory inference. Today we stop at the representation that makes these tasks possible. The examples come from mouse intestinal epithelium, so the same dataset will connect the stages. This is a review and tutorial, rather than a benchmark proving one universal pipeline. Its most useful distinction for today is between preparing a representation and drawing biological conclusions from it.',
    fa='شکل ۱ نقشهٔ کل مقاله است. امروز بخش چپ شکل و منطق رسیدن به UMAP را باز می‌کنیم. سمت راست یادآوری می‌کند که دیدن خوشه تازه شروع تفسیر زیستی است. این مقاله مرور آموزشی است و برندهٔ قطعی همهٔ روش‌ها را تعیین نمی‌کند.')

add('A cell begins as a vector of molecular counts','table',60,
    values=[['Cell','CD3D','MS4A1','NKG7','Total UMIs'],['A','4','0','1','1,000'],['B','8','0','2','2,000'],['C','0','6','0','1,200']],
    widths=[180,190,190,190,260],
    lines=['Each gene supplies one coordinate','Zero means no detected molecule','The full matrix contains many more genes'],
    note='Illustrative counts. Total UMIs include genes omitted from this table.',
    talk='Each row here is one cell and each gene is one coordinate. The table is intentionally tiny, but the actual matrix may contain thousands of cells and many thousands of genes. A count measures detected molecules in a UMI protocol, rather than a direct concentration. A zero can reflect lack of expression or limited sampling. Cells A and B have the same relative pattern in the displayed genes, but B has twice the library size. If we compare raw counts directly, sequencing depth can therefore influence which cells appear similar before we have learned anything about their identity.',
    fa='هر سلول یک بردار در فضای ژن‌هاست. مثال عمداً ساختگی است. A و B از نظر نسبت بیان ژن‌های نشان‌داده‌شده مشابه‌اند ولی عمق شمارش B دو برابر است. صفر الزاماً اثبات خاموش‌بودن ژن نیست. جهت ماتریس در نرم‌افزارها ممکن است برعکس باشد، اما مفهوم یکسان است.')

fig('QC changes which cells enter the geometry',2,75,
    points=['Count depth and detected genes','Mitochondrial fraction in context','Thresholds depend on tissue and sample'],
    talk='Figure 2 shows the first decision that changes our eventual map: which barcodes we retain. Panels A and B show count and gene distributions, and panel C shows the ranked count curve. Panel D combines detected genes, count depth and mitochondrial fraction. The red cutoffs belong to this particular intestinal dataset. They are examples, rather than universal thresholds. A high count can reflect a large cell or a doublet, and a high mitochondrial fraction needs tissue context. Removing a cell here removes its contribution to every later covariance, principal component and neighborhood relationship.',
    fa='در QC فقط داده را تمیز نمی‌کنیم؛ جمعیت سلول‌هایی را تعیین می‌کنیم که هندسهٔ بعدی بر اساس آن‌ها ساخته می‌شود. آستانه‌های قرمز شکل مخصوص همین داده‌اند. کم‌بودن RNA می‌تواند ویژگی یک جمعیت واقعی باشد، بنابراین تعداد ژن، UMI و درصد میتوکندری را با هم و در زمینهٔ هر نمونه می‌خوانیم.')

add('Normalization changes the meaning of a distance','table',75,
    values=[['Cell','Gene count','Total UMIs','Per 10,000','log1p'],['A','4','1,000','40','3.714'],['B','8','2,000','40','3.714']],
    lines=['Library-size normalization reduces depth differences','log1p compresses the expression range'],
    formula='yᵢⱼ = log(1 + 10,000 × cᵢⱼ / totalᵢ)',
    note='Illustrative arithmetic. Equal normalized values do not establish equal absolute RNA content.',
    source=PAPER+'\n'+SEURAT,
    talk='Library-size normalization divides a gene count by the total count in that cell, then applies a common scale factor. A and B now both have a normalized value of forty for this gene, and log one plus forty is approximately three point seven one four. The log transformation compresses the contribution of very large values. This changes the geometry on which PCA operates. The procedure assumes that relative expression is a useful comparison, and it does not recover absolute RNA content. Model-based alternatives, including SCTransform or residual-based workflows, change the normalization route but still require us to understand what representation enters dimensionality reduction.',
    fa='این مثال عددی نشان می‌دهد نرمال‌سازی چه چیزی را برابر می‌کند: سهم نسبی یک ژن در کتابخانه. لگاریتم دامنهٔ اعداد بزرگ را فشرده می‌کند. بیان نسبی با تعداد مطلق RNA برابر نیست. مسیر کلاسیک log-normalization را برای آموزش دنبال می‌کنیم؛ SCTransform یک مسیر جایگزین است.')

add('HVG selection chooses genes, scaling changes their weight','table',75,
    values=[['Operation','What changes','What stays meaningful'],['HVG selection','Columns used for geometry','Original counts remain available'],['Centering','Gene mean becomes zero','Relative position of cells'],['Unit-variance scaling','Each gene has variance one','Gene-wise standardized deviations']],
    widths=[240,400,460],
    lines=['HVGs emphasize variation beyond the mean–variance trend','Scaling is a modeling choice, not a universal requirement'],
    note='Classic log-normalized workflow. HVG methods differ in expected input.',
    talk='Feature selection and scaling do different jobs. HVG selection chooses which gene columns will define our representation, ideally emphasizing variability beyond the expected mean–variance relationship. Scaling changes the relative contribution of the selected genes. Centering removes a gene-specific mean, while unit-variance scaling expresses deviations in standard-deviation units. Scaling can prevent highly expressed genes from dominating, but it can also amplify noisy low-expression features. The original count matrix should remain available. The order and input depend on the HVG method, because some methods expect counts and others expect log-normalized expression. The choices here determine what variance PCA will regard as important.',
    fa='انتخاب HVG یعنی انتخاب محورهای ژنیِ ورودی؛ scaling یعنی تغییر وزن همان محورها. حذف ژن از محاسبهٔ هندسه نباید به معنی از دست‌دادن شمارش اصلی برای تحلیل بعدی باشد. روش‌های HVG ورودی یکسانی ندارند و ترتیب نسخهٔ کلاسیک را نباید به همهٔ روش‌ها تعمیم داد.')

fig('Batch correction can change the neighborhood structure',3,60,
    points=['Before and after correction','Mixing alone does not validate integration','Biology can be confounded with batch'],
    talk='Figure 3 compares a UMAP before and after batch correction, with cells colored by sample. The sample-specific arrangement changes after correction. This demonstrates that our map depends on how we define the representation before visualization. It does not prove that the corrected biology is right. We still need to check whether comparable cell identities align and whether real condition-specific states remain distinct. When condition and batch are completely confounded, a computational correction cannot establish which difference is technical. In current workflows, integration may adjust expression, a latent representation or a graph, so it does not have one mandatory position in every pipeline.',
    fa='رنگ‌ها نمونه‌ها را نشان می‌دهند. مخلوط‌شدن رنگ‌ها بعد از correction به‌تنهایی معیار موفقیت نیست؛ باید هویت‌ها و حالت‌های زیستی حفظ شوند. اگر گروه زیستی با batch کاملاً یکی باشد، الگوریتم نمی‌تواند منشأ اختلاف را اثبات کند. جای integration به نوع روش وابسته است.')

add('Many genes share a smaller number of variation patterns','text',60,
    lines=['Co-regulated genes contribute redundant information','Technical variation adds unwanted directions','A compact representation can improve neighbor search'],
    big='20,000 genes     2,000 HVGs     30 PCs     2 UMAP coordinates',
    note='Illustrative dimensions, not recommended fixed settings',
    talk='A cell has thousands of gene measurements, but many of those measurements vary together. A coordinated expression program can therefore occupy fewer effective dimensions than its number of genes. We also observe technical variation and measurement noise. PCA combines correlated features into a compact set of coordinates and can reduce the influence of some noise when we truncate the representation. The numbers on this slide describe an illustrative route, not fixed settings. Retaining thirty PCs is still a thirty-dimensional analysis. The final move to two UMAP coordinates serves the separate goal of making the relationships visible on a screen.',
    fa='تعداد ژن با تعداد الگوهای مستقل تغییرات یکی نیست. یک برنامهٔ زیستی ممکن است ده‌ها ژن را هم‌جهت تغییر دهد. PCA از هم‌بستگی‌ها برای خلاصه‌کردن استفاده می‌کند. عدد ۳۰ فقط مثال است؛ قبل از UMAP هنوز در یک فضای ۳۰بعدی مقایسه می‌کنیم، نه روی دو محور اول PCA.')

add('PCA rotates the axes, truncation removes directions','wideimage',90,
    image=str(A/'pca-geometry.png'),
    caption='PC1 captures the largest remaining variance. PC2 is orthogonal to PC1.',
    source=PAPER+'\nConceptual educational image generated with OpenAI ImageGen. Illustrative geometry, no experimental measurements.',
    talk='In the left plot, the two genes vary together, so most of the spread lies along a diagonal direction. PCA defines PC1 along the direction with the greatest variance and PC2 along the perpendicular direction with the next greatest variance. The right plot expresses the cloud in these new coordinates. Keeping both components is a rotation and preserves the full information in this two-dimensional example. Dimensionality reduction occurs when we retain PC1 and discard the smaller perpendicular variation. The orange ticks illustrate that projection. In real data, discarded variation can include noise and useful biology, which is why truncation needs inspection.',
    fa='خود چرخشِ همهٔ محورها چیزی را حذف نمی‌کند. کاهش بُعد وقتی اتفاق می‌افتد که فقط چند PC را نگه داریم. PC1 بیشترین واریانس را می‌گیرد و PC2 در جهت عمود، بیشترین واریانس باقیمانده را. مقدار کمِ واریانس معادل بی‌ارزش‌بودن زیستی نیست؛ یک جمعیت نادر ممکن است در مؤلفه‌های بعدی دیده شود.')

add('Loadings describe genes, scores locate cells','text',90,
    formula='T = ZVₖ',
    lines=['Z: centered expression for the selected genes','Vₖ: gene weights for each retained component','T: one row per cell, one column per PC'],
    big='PC1 score = 0.6 × Gene A + 0.8 × Gene B',
    note='Illustrative loading vector. Genes A and B denote centered, optionally scaled expression.',
    talk='The matrix equation separates two objects that are often confused. Z contains the centered, and possibly scaled, expression values. The columns of V are loading vectors, with one weight per gene. Multiplying them gives T, the scores that locate each cell in PC space. In the simple example, a cell high in both genes receives a high PC1 score. A loading can also be negative, so a component may contrast two expression programs. A principal component is therefore a weighted combination of genes rather than one selected gene or an assigned cell type. Loading inspection connects the numerical direction back to possible biology.',
    fa='loading وزن ژن در یک مؤلفه است و score مختصات یک سلول روی همان مؤلفه. از روی ژن‌های با وزن مثبت و منفی می‌توان معنای احتمالی PC را بررسی کرد. PC معادل یک ژن یا یک خوشه نیست. تغییر علامت کل یک PC هیچ معنای زیستی تازه‌ای ایجاد نمی‌کند.')

add('Truncating PCA changes the distances between cells','text',75,
    formula='d²(i,j) = Σₗ₌₁ᵏ (Tᵢₗ − Tⱼₗ)²',
    lines=['The retained PCs define the comparison space','All PCs preserve centered Euclidean geometry','Discarded PCs no longer contribute to distance'],
    note='For an orthonormal PCA basis, without whitening',
    talk='The retained scores define the distances used for neighbor search. With an orthonormal PCA basis and every available component retained, PCA preserves Euclidean distances in the centered input space. Once we truncate, only differences along the retained components contribute. This explains both the benefit and the risk. If discarded directions mainly contain noise, neighborhoods can become more useful. If a rare cell state differs along a discarded direction, its separation can disappear. Standard PCA scores retain their variance differences. Whitening would reweight the components and change the geometry again, so it should not be confused with ordinary PCA.',
    fa='فاصلهٔ دو سلول در فضای PCA از اختلاف scoreها در مؤلفه‌های نگه‌داشته‌شده محاسبه می‌شود. با نگه‌داشتن همهٔ مؤلفه‌ها، چرخش فاصله را عوض نمی‌کند؛ با حذف مؤلفه‌ها، سهم آن جهت‌ها حذف می‌شود. بنابراین PCA تضمین پاک‌کردن نویز نیست و ممکن است بخشی از سیگنال را هم حذف کند.')

add('The elbow is a starting point for choosing PCs','figure',75,
    image=str(A/'Fig4F.png'),label='Figure 4F',legend=legend(paper,'Figure 4'),
    points=['This dataset has an elbow around PCs 5–7','Inspect loadings and QC associations','Check whether neighborhoods remain stable'],
    talk='Panel F of Figure 4 ranks the variance explained by the first thirty-one principal components. The steep decline becomes flatter around PCs five to seven in this intestinal dataset. The paper uses that elbow as a practical heuristic. It does not establish a universal PC count. A better decision also considers the genes loading on the components, associations with depth or batch, and the stability of relevant neighborhoods over plausible choices. A low-variance component may describe a small population. Conversely, a dominant component may reflect technical variation. The goal is a useful representation for the biological question, rather than maximizing retained variance alone.',
    fa='آرنج نمودار برای این داده حدود PCهای ۵ تا ۷ است. این عدد را برای دادهٔ خودمان کپی نمی‌کنیم. loadings، رابطه با QC و batch، و پایداری همسایگی‌ها را در چند انتخاب معقول بررسی می‌کنیم. بیشترین واریانس الزاماً مهم‌ترین زیست‌شناسی نیست.')

add('A strong PC can reflect sequencing depth','figure',60,
    image=str(A/'Fig4A.png'),label='Figure 4A',legend=legend(paper,'Figure 4'),
    points=['Cells share the count-depth color scale in Figure 4','PCA optimizes variance without biological labels','Depth associations require investigation'],
    talk='Panel A shows the first two PCs of the same intestinal cells. The original figure colors cells by count depth, and high-depth cells occupy a distinct part of the representation. PCA has no knowledge of cell types or technical quality. It simply finds variation in the matrix we supplied. A depth association can indicate residual technical structure, but depth can also correlate with real cell size or activity. This is a reason to inspect the relationship and the relevant genes, rather than automatically delete PC1. The question is whether the component captures the biology we want to compare.',
    fa='در شکل اصلی، رنگ نقاط عمق شمارش است. PCA برچسب زیستی نمی‌شناسد و حتی PC1 می‌تواند با عمق یا batch همراه باشد. این مشاهده دلیل بررسی است، نه مجوز حذف خودکار PC1. عمق گاهی با اندازه یا فعالیت واقعی سلول مرتبط است.')

add('Neighbor search turns PC coordinates into local relationships','table',75,
    values=[['Candidate','PC coordinates','Distance from A','Among 2 nearest?'],['A','(0, 0, 0)','0','Reference'],['B','(0.1, 0.2, 0.1)','0.245','Yes'],['C','(2, 0, 0)','2','Yes'],['D','(0, 3, 0)','3','No']],
    widths=[170,390,260,290],
    lines=['Each cell searches in the retained PC space','A graph records local connectivity and weights'],
    note='Toy example with Euclidean distance and k = 2. A is excluded from its own candidates.',
    source=PAPER+'\n'+ALGO,
    talk='The next object is a neighborhood graph. In this toy example, A is the reference cell and B and C are its two nearest neighbors in three-dimensional PC space. In a real analysis, the same search uses all retained PCs, not just the first two on a plot. Repeating the search for each cell produces local relationships. UMAP turns local distances into weighted, fuzzy connections, adapting to local scales. This is why the number of neighbors matters. A graph emphasizes who is connected to whom, and does not retain every pairwise distance with equal importance.',
    fa='گام واسطِ فراموش‌شده بین PCA و UMAP ساخت رابطهٔ همسایگی است. جدول مثال ساختگی سه‌بعدی است و سلول A همسایهٔ خودش حساب نمی‌شود. در دادهٔ واقعی از همهٔ PCهای انتخاب‌شده استفاده می‌شود. وزن یال در UMAP صرفاً خود فاصلهٔ خام نیست و به مقیاس محلی هم بستگی دارد.')

add('UMAP fits a two-dimensional layout to a weighted graph','text',90,
    lines=['Local distances define connection strengths','Nearby cells attract in the embedding objective','Repulsion helps separate weakly connected cells','Optimization trades off conflicting relationships'],
    big='The output gives each cell two display coordinates',
    note='Conceptual summary of UMAP’s fuzzy graph objective and negative-sampling optimization',
    source=ALGO,
    talk='UMAP first builds a weighted description of local relationships in the input space. It then searches for low-dimensional coordinates whose connection strengths resemble that description. Strong connections create attraction in the optimization, while negative sampling approximates a repulsive contribution for other pairs. A flat map cannot satisfy all relationships in a complex dataset, so the algorithm makes tradeoffs. Initialization, the optimization and the input graph all affect the arrangement. The axes are coordinates of this fitted layout, and do not have PCA-like gene loadings. UMAP can embed into more dimensions, but here we use two because our purpose is visualization.',
    fa='UMAP ابتدا رابطه‌های محلی وزن‌دار می‌سازد و سپس مختصاتی می‌یابد که آن رابطه‌ها را در بُعد کم تقریب بزند. نمایش تخت ناچار به مصالحه است. محور UMAP1 وزن ژنی قابل‌خواندن مثل یک PC ندارد. UMAP الزاماً محدود به دو بُعد نیست، ولی اینجا دو بُعد را برای نمایش می‌خواهیم.')

fig('The same cells produce different views of structure',4,75,
    wide=True,
    talk='Figure 4 compares several views of the same cells, with count depth supplying the shared color scale. Panel A uses the first two PCs. Panels B through E show t-SNE, diffusion maps, UMAP and a force-directed layout. Panel F gives the PCA variance curve. The shapes differ because the methods optimize different properties. A compact PCA representation can still contain useful structure even when its first two coordinates overlap visually. UMAP can make local relationships easier to see in two dimensions, but the appearance does not independently validate the cell states or guarantee faithful distances between separated groups.',
    fa='همان سلول‌ها و همان متغیر رنگ در چند روش، شکل‌های متفاوتی می‌سازند. نامرتب‌بودن دو PC اول به معنی بی‌فایده‌بودن فضای چندده‌بعدی PCA نیست. UMAP برای نمایش دوبعدی مناسب‌تر است، اما تصویر جذاب‌تر صحت زیستی یا فاصلهٔ دقیق بین جزیره‌ها را ثابت نمی‌کند.')

add('min_dist changes how tightly points pack','pairimage',60,
    images=[str(A/'min-dist-0.png'),str(A/'min-dist-099.png')],
    caption='Same synthetic color dataset from the official UMAP tutorial',
    lines=['Smaller min_dist allows tighter packing','This changes the layout, not the input expression'],
    source=UMAP,
    talk='These two official tutorial plots use the same synthetic color dataset and change the minimum-distance setting. The smaller value permits more tightly packed groups, while the larger value spreads points more broadly. The example is deliberately not a cell dataset, which helps separate an algorithmic effect from a biological interpretation. In a workflow with a fixed input graph, changing min_dist changes the display rather than the expression matrix. The number of neighbors has a different role because it changes the local relationships considered. A more compact island should therefore prompt a parameter check before a claim about a sharper cell identity.',
    fa='هر دو تصویر از یک دادهٔ مصنوعی رنگ‌ها هستند، نه سلول واقعی. min_dist کوچک‌تر اجازهٔ تجمع فشرده‌تر می‌دهد. با گراف ثابت، این تنظیم ظاهر نقشه را تغییر می‌دهد. n_neighbors نقش دیگری دارد و تعریف همسایگی را تغییر می‌دهد؛ این دو پارامتر قابل‌جایگزینی نیستند.')

add('PCA and UMAP solve different problems','table',75,
    values=[['Property','PCA used for representation','UMAP used for display'],['Objective','Capture variance in orthogonal axes','Approximate weighted local relationships'],['Typical output here','Several dozen coordinates','Two coordinates'],['Interpretation','Gene loadings and cell scores','Relative layout and local neighborhoods'],['Main limitation','Linear axes may miss nonlinear structure','Global distances and density can distort']],
    widths=[220,440,460],
    note='UMAP can run without PCA and can output more than two dimensions.',
    talk='The sequence makes sense because the objectives differ. PCA provides a compact linear representation with interpretable gene loadings, and UMAP makes a local relationship structure visible in two dimensions. Applying UMAP directly to all genes is possible, so PCA is not a mathematical prerequisite. In typical single-cell workflows, PCA offers computational efficiency and a way to reduce redundant or noisy directions before neighbor search. These benefits depend on the preprocessing and the retained components. UMAP cannot restore a rare-state signal that we removed earlier. Two reductions are useful when their different objectives match our analysis and communication needs.',
    fa='PCA پیش‌شرط ریاضی UMAP نیست. استفادهٔ رایج از این زنجیره به دلیل تفاوت هدف‌هاست: ساخت یک نمایش فشرده و بعد قابل‌دیدن‌کردن روابط آن. اگر در HVG یا حذف PCها سیگنالی را دور ریخته باشیم، UMAP آن را برنمی‌گرداند.')

add('Clustering and UMAP branch from related representations','table',60,
    values=[['Task','Usual input','Result'],['Neighbor construction','Retained PCs or integrated latent space','Weighted cell relationships'],['Clustering','A neighborhood graph, often with SNN weighting','A discrete label per cell'],['UMAP visualization','PC-derived neighborhoods or a supplied graph','Two display coordinates']],
    widths=[280,520,310],
    lines=['Cluster labels usually appear as colors on UMAP','Changing layout settings need not change the clusters'],
    source=PAPER+'\n'+SCAN+'\n'+SEURAT,
    talk='A common misunderstanding is that the visible UMAP islands are the input to routine single-cell clustering. In typical workflows, graph clustering and visualization branch from related high-dimensional representations. Scanpy can use a neighborhood graph for both Leiden clustering and UMAP. Seurat commonly uses an SNN graph for clustering while RunUMAP builds its own neighborhoods from the specified PCs. These graphs need not be identical. Once labels have been assigned, we color the UMAP by those labels. Changing the layout alone can therefore alter the picture without changing the labels. Other clustering strategies exist, but this is the standard distinction relevant here.',
    fa='معمولاً خوشه‌بندی روی جزیره‌های دوبعدی UMAP انجام نمی‌شود. در Scanpy گراف همسایگی می‌تواند برای هر دو استفاده شود. در Seurat گراف SNN خوشه‌بندی و همسایگی داخلی RunUMAP الزاماً یکی نیستند. رنگ‌کردن نقاط با برچسب خوشه یعنی نتیجهٔ تحلیل دیگری را روی نقشه نشان می‌دهیم.')

add('UMAP distances need evidence beyond the picture','text',60,
    lines=['Nearby cells suggest local similarity worth checking','A gap does not prove a discrete cell type','Island area does not estimate population abundance','An apparent path does not establish lineage or time'],
    big='Would the conclusion survive a different reasonable representation?',
    talk='The map is useful for generating questions. Nearby cells may share local expression structure, but we should check the relevant markers and representation. A gap can depend on preprocessing, neighborhood choices and the layout. The area of an island does not estimate the fraction of cells in that population, because UMAP does not preserve density in a straightforward way. An apparent bridge does not establish direction, developmental time or lineage. Strong biological claims need independent evidence and a design that supports them. A useful journal-club question is whether the conclusion survives plausible choices of features, PCs and neighborhood settings.',
    fa='از شکل جزیره نباید فراوانی جمعیت، مرز قطعی نوع سلولی یا مسیر زمانی استنتاج کنیم. اندازهٔ جزیره با تعداد نقاط یکی نیست و UMAP چگالی را ساده و مستقیم حفظ نمی‌کند. نتیجهٔ زیستی را با مارکرها، طراحی مطالعه، تکرارهای زیستی و پایداری نمایش بررسی می‌کنیم.')

add('Different analysis tasks use different data layers','table',45,
    values=[['Task','Appropriate representation'],['Library-aware count models and pseudobulk','Retained raw counts with sample metadata'],['Gene-expression display','An appropriate normalized gene-expression layer'],['Neighbor search and graph clustering','Selected PC scores or a suitable latent space'],['Visual exploration','UMAP coordinates, colored by independent metadata']],
    widths=[590,520],
    note='Teaching summary based on Table 1. Statistical input depends on the chosen method.',
    legend=legend(paper,'Table 1'),label='Table 1, teaching summary',
    talk='Table 1 in the paper reminds us that analysis produces multiple data layers rather than a single progressively better matrix. We preserve raw counts because library-aware models and sample-level pseudobulk analyses need them. We use appropriate gene-expression values to display markers, and PC scores or another latent representation to compare cells. UMAP coordinates are mainly a visualization product in this workflow. Differential-expression methods have different input requirements, and a biological comparison also requires the correct experimental unit. The practical lesson is to retain provenance and choose the representation for the task instead of feeding the latest output into every subsequent analysis.',
    fa='این جدول برداشت آموزشی از Table 1 است، نه بازنویسی لفظی همهٔ توصیه‌های سال ۲۰۱۹. یک ماتریسِ نهایی برای همهٔ کارها نداریم. شمارش خام، بیان نرمال‌شده، PCA و UMAP هر کدام کاربرد خود را دارند. برای مقایسهٔ شرایط، واحد آزمایشی و تکرار زیستی هم ضروری‌اند.')

add('The full chain explains the final map','text',45,
    lines=['Counts record sampled molecules','QC and normalization define comparable observations','HVGs and scaling define the gene space','PCA and neighbor search define local similarity','UMAP supplies the display coordinates'],
    big='Every visible island inherits upstream decisions',
    talk='We can now answer the title question. Counts record sampled molecules. Quality control chooses the cells we analyze, and normalization defines a useful comparison scale. Feature selection and scaling decide which gene differences matter. PCA summarizes that representation, and neighbor search expresses local similarity. UMAP fits display coordinates to those relationships. This order is a useful classical route, with method-specific alternatives at several steps. The final map inherits each upstream decision. Understanding the sequence allows us to diagnose a surprising plot by returning to the representation that produced it, rather than adjusting the appearance alone.',
    fa='زنجیره را یک بار با زبان ساده جمع می‌کنیم. هر مرحله نوعی تصمیم دربارهٔ مشاهده، مقایسه یا نمایش است. این ترتیب مسیر کلاسیک آموزشی است و بسته به روش می‌تواند تغییر کند. برای عیب‌یابی یک UMAP عجیب باید به ورودی و همسایگی برگردیم، نه فقط تنظیم ظاهر.')

add('Discussion: which decisions would change your conclusion?','discussion',300,
    lines=['A rare population disappears when we use fewer PCs','A cluster becomes tighter after changing min_dist','Samples mix after integration, but a disease state disappears'],
    note='Which result needs a representation check, a layout check, or biological validation?',
    talk='The paper gives us a strong organizing framework, especially the distinction between summarization and visualization. Its limitations are equally useful for discussion. It is a 2019 review, its worked example is one tissue system, and a heuristic such as an elbow does not guarantee biological completeness. Consider these three observations. A rare population vanishes when we retain fewer PCs. A cluster becomes tighter after changing min_dist. Integration mixes samples but removes a disease-associated state. In each case, which object changed, what evidence would we inspect next, and what conclusion can we support without treating the UMAP picture as ground truth?',
    fa='برای پنج دقیقهٔ آخر، سه سناریو را مقایسه می‌کنیم. ناپدیدشدن جمعیت با کاهش PC مسئلهٔ نمایشِ تحلیلی است؛ فشرده‌شدن با min_dist بیشتر مسئلهٔ چیدمان است؛ حذف حالت بیماری بعد از integration نیازمند ارزیابی حفظ زیست‌شناسی است. مقاله چارچوب فکری مفیدی دارد، ولی مرور ۲۰۱۹ و یک مثال بافتی تضمین نسخهٔ عمومی نیستند.')

add('Appendix: complete figures and implementation','divider',0,
    subtitle='Additional evidence for questions after the 30-minute session',
    talk='The remaining slides retain the downstream figures, expanded-view evidence, the appendix example and implementation references for discussion.',
    fa='از اینجا اسلایدهای پشتیبان شروع می‌شوند و در زمان اصلی ارائه قرار ندارند. تمام شکل‌های باقی‌ماندهٔ مقاله و ضمائم و دو مثال کد در این بخش‌اند.')

fig('Downstream tasks ask different biological questions',5,0,wide=True,
    talk='Figure 5 groups downstream tasks by the level at which they ask a question. The cell-level branch includes clustering and trajectory analysis, while the gene-level branch includes differential expression, gene sets and regulatory networks. Cluster analysis also includes annotation and comparisons of composition. The same UMAP-like arrangement appears under several tasks, but the picture is only a common display. The appropriate input, statistical unit and interpretation differ across tasks. This figure explains why reaching a visually organized map is an intermediate result. We still need a defined biological question and a method whose assumptions match that question.',
    fa='شکل ۵ کاربردهای پایین‌دستی را در سطح سلول و ژن جدا می‌کند. استفاده از یک نقشهٔ مشابه برای چند سؤال به معنی یکسان‌بودن ورودی یا روش آماری آن‌ها نیست. خوشه، trajectory، ترکیب جمعیتی و بیان تفاضلی هر کدام سؤال متفاوتی دارند.')

fig('Markers support the interpretation of cluster labels',6,0,wide=True,
    talk='Figure 6 connects graph-derived clusters to biological interpretation. Panel A shows annotated cell identities in a UMAP representation. Panel B overlays example marker genes, including Slc12a2 for stem cells, Arg2 for enterocytes, Tff3 for goblet cells and Defa24 for Paneth cells. Panel C compares the distribution of proximal and distal enterocytes in the map. The clustering supplies groups, while marker expression and prior biological knowledge help interpret those groups. A separated island alone cannot establish a cell type. Expression overlays should come from an appropriate gene-expression layer, and the labels need supporting evidence beyond the visual separation.',
    fa='شکل ۶ از رنگ خوشه به معنای زیستی آن می‌رسد. مارکرها برای تفسیر برچسب‌ها کمک می‌کنند، نه اینکه جزیرهٔ UMAP به‌تنهایی نوع سلولی را اثبات کند. پنل C محل enterocyteهای proximal و distal را مقایسه می‌کند. پنل‌ها باید همراه با legend و شناخت بافت خوانده شوند.')

fig('Trajectory inference adds assumptions beyond a UMAP layout',7,0,wide=True,
    talk='Figure 7 shows several representations used to study a continuous process. Panel A maps enterocyte relationships, and panel B fits Slingshot trajectories through clusters in PCA space. Panel C relates cluster distributions to pseudotime. Panel D shows a graph abstraction projected onto a UMAP, and panel E displays gene-expression patterns over the inferred ordering. These panels combine a representation with additional modeling assumptions. The trajectory does not emerge merely because a UMAP contains an apparent bridge. Root selection, connectivity and the biological interpretation of pseudotime all matter. Cross-sectional transcriptional similarity alone does not directly measure elapsed time or prove lineage.',
    fa='شکل ۷ نشان می‌دهد مسیر و pseudotime از مدل و فرض‌های اضافه می‌آیند. Slingshot در پنل B در فضای PCA کار می‌کند و پنل D گراف انتزاع‌شده را روی UMAP نمایش می‌دهد. شباهت رونویسی یا یک پل ظاهری، به‌تنهایی زمان واقعی یا دودمان را اثبات نمی‌کند.')

for n,title,talk,fa in [
    (1,'Published workflows used a range of HVG counts','Figure EV1 summarizes a manual survey of published workflows. The horizontal axis gives the logarithm of cell number, the vertical axis gives the number of highly variable genes, and colors indicate technology. The survey shows considerable variation in choices across studies. It provides context for the paper’s suggested range, rather than a controlled test of the optimal number of genes. Dataset size alone does not identify a universally correct setting. Tissue complexity, the signal of interest and the feature-selection method also matter. The underlying survey records are supplied in Dataset EV1.','این شکل پیمایش دستی مقاله‌هاست و benchmark کنترل‌شده نیست. تنوع انتخاب تعداد HVG را نشان می‌دهد. از روی تعداد سلول‌ها به‌تنهایی نمی‌توان بهترین تعداد ژن را تعیین کرد.'),
    (2,'Correction and denoising change gene-level variation','Figure EV2 shows changes in the coefficient of variation after correction. Panels A and B use ComBat batch correction for intestinal epithelium and embryonic stem-cell data. Panels C and D use DCA denoising for the corresponding datasets. Negative changes indicate a decrease in the coefficient of variation. The distributions make clear that processing can alter gene-level variation substantially. They do not by themselves establish which changes remove technical noise and which affect useful biology. This matters for our central topic because variance-changing operations upstream can alter the principal components and the neighborhoods used for visualization.','ردیف بالا تغییر CoV بعد از ComBat و ردیف پایین بعد از DCA را نشان می‌دهد. مقادیر منفی یعنی کاهش CoV. تغییر واریانس به‌تنهایی به معنی بهبود نیست؛ باید معلوم شود نویز فنی کم شده یا سیگنال زیستی آسیب دیده است.')]:
    add(title,'figure',0,image=str(A/f'FigEV{n}.png'),label=f'Figure EV{n}',legend=legend(ev,f'Figure EV{n}'),wide=True,talk=talk,fa=fa)

add('Clustering random data can create apparent marker evidence','figure',0,
    image=str(A/'AppendixFig1.png'),label='Appendix Figure 1',wide=True,
    legend='Appendix Figure 1: Analysis of marker gene detection in clustered random data. Random data were simulated using the splatter package for one cluster, with dropout, and without differentially expressed genes. After Louvain community detection we detected 2 clusters, which are visualized in PCA-space in a. Marker gene detection was performed in scanpy using the t-test. The distribution of the p-values of all genes is shown in b.',
    talk='The appendix provides a useful warning about circular validation. The authors simulated a single population without differentially expressed genes, yet graph clustering produced two groups. Panel A shows those groups in PCA space. Panel B shows an excess of small p-values when marker testing uses the same data that helped define the groups. The test ignores the selection created by clustering. A visually distinct arrangement and significant markers from the same data therefore cannot serve as independent proof of a biological population. External validation or a design that accounts for the selection is needed for stronger claims.',
    fa='دادهٔ شبیه‌سازی‌شده یک جمعیت دارد و DE واقعی ندارد، اما خوشه‌بندی دو گروه پیدا می‌کند. آزمون مارکر روی همان داده می‌تواند p-value کوچک بدهد. این اعتبارسنجی مستقل نیست و انتخاب ناشی از خوشه‌بندی را نادیده می‌گیرد. نکته برای نقد ادعاهای زیستی از روی نقشه مهم است.')

add('The paper assigns different uses to each processing stage','table',0,
    values=[['Layer','Stage','Applications in the 2019 paper'],['Measured','1 Raw\n2 Normalized (+ log transformed)','Statistical testing'],['Corrected','3.1 Technical correction','Visual comparison'],['Corrected','3.2 Biological correction','Preprocessing for trajectory inference'],['Reduced','4 Feature selected','Visualization, trajectory inference'],['Reduced','5 Dimensionally reduced (summarized)','Visualization, clustering, KNN graphs, trajectory inference']],
    widths=[190,390,530],label='Table 1',legend=legend(paper,'Table 1'),
    note='Full stage coverage, wording condensed. Method-specific assumptions still apply.',
    talk='This editable table retains the stages and application groups listed in Table 1 of the original paper. Measured, corrected and reduced layers serve different purposes. The authors distinguish technical correction from removal of biological covariates, and distinguish feature selection from dimensional summarization. These categories are useful, but should be interpreted in the context of the specific method. A current count model may require raw counts, and a modern integration model may produce a latent representation rather than corrected expression. The enduring lesson is to select the data layer deliberately and retain the information required by each downstream analysis.',
    fa='در این جدول پوشش کامل مراحل Table 1 مقاله حفظ شده و متن کوتاه شده است. تقسیم‌بندی measured، corrected و reduced مفید است، ولی کاربرد دقیق به فرض روش وابسته است. پیشنهادهای ۲۰۱۹ را به همهٔ مدل‌های جدید تعمیم نمی‌دهیم.')

add('A classic Seurat route makes the representations explicit','code',0,
    code='obj <- NormalizeData(obj)\nobj <- FindVariableFeatures(obj, nfeatures = 2000)\nobj <- ScaleData(obj, features = VariableFeatures(obj))\nobj <- RunPCA(obj, features = VariableFeatures(obj))\n\nobj <- FindNeighbors(obj, dims = 1:30)\nobj <- FindClusters(obj, resolution = 0.5)\nobj <- RunUMAP(obj, reduction = "pca", dims = 1:30)',
    note='Assumes QC is complete. Illustrative settings, not optimized defaults. SCTransform follows an alternative route.',
    source=SEURAT,
    talk='This code spells out the classical log-normalized route in Seurat. Variable features define the genes entering PCA, and the explicit dimensions identify the PC coordinates used downstream. FindNeighbors constructs the graph used for clustering. RunUMAP, when given a reduction and dimensions this way, computes its own neighborhood-based embedding from those coordinates. The final labels can then be displayed on the map. The numbers here are illustrative, and the same chosen dimensions should be justified for the dataset. This snippet assumes quality control is complete and preserves the original counts in the object. SCTransform provides an alternative preprocessing route.',
    fa='کد برای توضیح اتصال مراحل است و روی دادهٔ کاربر اجرا نشده است. QC از قبل فرض شده. ۲۰۰۰ ژن، ۳۰ PC و resolution=0.5 نسخهٔ قطعی نیستند. RunUMAP در این شکل، مستقیماً گراف SNN خوشه‌بندی را دریافت نمی‌کند و از PCA همسایگی می‌سازد.')

add('A classic Scanpy route exposes the shared neighbor graph','code',0,
    code='adata.layers["counts"] = adata.X.copy()\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\nsc.pp.highly_variable_genes(adata, flavor="seurat",\n                          n_top_genes=2000)\nsc.pp.scale(adata, max_value=10)\nsc.pp.pca(adata, n_comps=50, mask_var="highly_variable")\nsc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)\nsc.tl.umap(adata, random_state=0)\nsc.tl.leiden(adata, resolution=0.5)',
    note='Assumes QC is complete and X initially contains counts. Illustrative settings. Scaling can increase memory use.',
    source=SCAN+'\nhttps://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.pca.html',
    talk='This Scanpy example keeps a count layer before normalization, then selects variable genes using a flavor that expects log-normalized expression. PCA uses the explicit highly-variable mask. The neighbor step creates the graph that UMAP and Leiden can both use, which makes the branch visible in the code. Running the visualization before or after assigning cluster labels does not make those labels depend on the two-dimensional layout. Scaling can increase memory use, so this compact teaching example should be adapted for large sparse data. These are illustrative settings rather than a fitted analysis, and the dependencies and input assumptions still require checking.',
    fa='برای flavor="seurat" ورودی log-normalized است؛ flavor="seurat_v3" انتظار شمارش دارد و جایش فرق می‌کند. گراف neighbors می‌تواند ورودی مشترک UMAP و Leiden باشد. نگه‌داشتن counts برای تحلیل بعدی مهم است. scaling ممکن است حافظهٔ زیادی مصرف کند و این کد مثال آموزشی است.')

add('Sources and the scope of the update','sources',0,
    lines=['Main paper: Luecken & Theis, MSB 2019, e8746','Main Figures 1–7, Expanded View Figures 1–2, Appendix Figure 1','Current implementation references: Seurat, Scanpy, UMAP','Added arithmetic and PCA geometry are teaching examples'],
    small='doi.org/10.15252/msb.20188746\nsatijalab.org/seurat/articles/pbmc3k_tutorial\nscanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html\numap-learn.readthedocs.io/en/latest/how_umap_works.html\numap-learn.readthedocs.io/en/latest/parameters.html',
    source=PAPER+'\n'+SEURAT+'\n'+SCAN+'\n'+ALGO+'\n'+UMAP,
    talk='The 2019 tutorial remains the scientific anchor for this presentation. All seven main figures, both expanded-view figures and the appendix figure are retained, with their source legends in the notes. The current documentation clarifies implementation details such as the neighborhood graph and the inputs expected by specific functions. The arithmetic tables and PCA geometry are illustrative teaching additions rather than new experimental results. Recommendations that depend on software or data type should be checked against the documentation and the biological question, while the distinction between summarization, neighborhood construction and visualization remains the central organizing idea.',
    fa='مرز منابع روشن است: شکل‌ها و چارچوب مفهومی از مقاله‌اند، توضیح کد از مستندات رسمی فعلی و مثال‌های عددی و هندسی ساختهٔ آموزشی‌اند. این ارائه benchmark تازه‌ای انجام نمی‌دهد و نتیجهٔ آزمایش جدیدی ادعا نمی‌کند.')

# Keep a precise 25-minute main route, followed by five minutes of discussion.
main=S[:22]
delta=1500-sum(x['seconds'] for x in main)
main[-1]['seconds']+=delta
assert sum(x['seconds'] for x in S[:23])==1800
for i,s in enumerate(S,1):
    s['number']=i
    s['appendix']=i>23
    s['notes']='Ready-to-read note:\n'+s['talk']
    if s.get('legend'):s['notes']+='\n\nLegend:\n'+s['legend']
    s['notes']+='\n\nSource:\n'+s['source']+'\n\nPersian explanation:\n'+s['fa']
    s['notes']+='\n\nTiming:\n'+('Appendix, outside the 30-minute route.' if s['appendix'] else f"{s['seconds']} seconds. Cumulative {sum(x['seconds'] for x in S[:i])//60:02}:{sum(x['seconds'] for x in S[:i])%60:02}.")
(ROOT/'build/deck-content.json').write_text(json.dumps(S,ensure_ascii=False,indent=2))
guide=['# From Count Matrix to UMAP: What Happens in Between?','',
       'Main route: slides 1–22, 25 minutes. Discussion: slide 23, 5 minutes. Slides 24 onward are backup material.','',
       'Scientific anchor: Luecken & Theis (2019), doi:10.15252/msb.20188746.','',
       'Visible slide text is English. Every slide has an English script and Persian preparation notes.','']
for s in S:
    guide += [f"## {s['number']}. {s['title'].replace(chr(10),' ')}",'',('Appendix' if s['appendix'] else f"{s['seconds']} seconds"),'',s['talk'],'',s['fa'],'']
(ROOT/'Presenter-guide.md').write_text('\n'.join(guide))
print('Slides',len(S),'Main timing',sum(s['seconds'] for s in S[:23]),'seconds','adjustment',delta)

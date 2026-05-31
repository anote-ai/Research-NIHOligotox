U.S. flag
An official website of the United States government

Here’s how you know
Home
Health Information
Grants & Funding
News & Events
Research & Training
Institutes at NIH
About NIH
Virtual Tour
En Español
Search
Search NIH

Search
Breadcrumb
Home  NIH Challenges and Prize Competitions
On this page
Overview
Timeline
Prizes
Rules
Judging
How to enter
Contact
Winners
Oligonucleotide Toxicity (OligoTox) Open Data Challenge
Generation and public dissemination of data on oligonucleotide toxicity in human cells.

Image
Oligonucleotide Toxicity (OligoTox) Open Data Challenge
The Oligonucleotide Toxicity Open Data Challenge incentivizes publicly accessible high-quality datasets from human cells to predict oligonucleotide toxicity from sequence and chemical modification.

Phase 2 open until 12/31/2026 11:59 PM ET

Total cash prizes: $500,000

 

Overview
Background

Oligonucleotides (“oligos”) are an emerging class of therapeutics most commonly used in the treatment of rare genetic diseases. As used here, the term “oligos” refers to short (typically 15-30 nucleotides in length) nucleic acids intended to be used as therapeutics to treat human diseases by modulating gene expression. Typically, oligo therapeutics contain chemical modifications in the phosphodiester backbone or sugar moiety to limit degradation by cellular nucleases. Examples of therapeutic oligos include anti-sense oligos (ASOs), siRNAs, and microRNA mimics. To date, there are 17 FDA-approved oligo therapeutics, and many more in clinical development, to treat rare genetic diseases, infectious diseases, and cancer.

As with any new molecular entity used as a therapeutic, the consideration of toxicity is a critical requirement in assessing safety. Toxicity in oligo therapeutics may be sequence-dependent (e.g., hybridization to nucleic acids other than the intended target) or sequence-independent (e.g., binding to cellular proteins) and are primarily assessed through expensive and time-consuming animal studies. Notably, both the FDA (FDA Announces Plan to Phase Out Animal Testing Requirement for Monoclonal Antibodies and Other Drugs) and NIH (NIH to prioritize human-based research technologies | National Institutes of Health) are focusing on new approach methodologies that can reduce animal testing in preclinical safety studies through practices such as use of in vitro human-based systems, in silico/computational modeling, and other innovative platforms that can collectively evaluate immunogenicity, toxicity, and pharmacodynamics in humans and improve the predictive relevance of preclinical drug testing. One possible strategy for assessing the toxicity of oligo therapeutics would be the use of in vitro human-based systems, including microphysiological systems (i.e., tissue chips, organs-on-a-chip) and 3-dimensional human organoids. NIH has already invested in both of these areas, and research is ongoing.

In recent years, we have witnessed dramatic increases in the use of advanced computational methods to address complex biological problems, such as AlphaFold for the accurate prediction of protein structures. Advances in data science, such as through machine learning (ML) and artificial intelligence (AI), and in computing ecosystems, such as with quantum computing, present an opportunity to explore new, innovative approaches that may lead to the development of safe and effective oligo therapeutics. An RFI on Better Tools to Predict Toxicities Resulting from Oligonucleotide Therapeutics was published by NCATS to solicit input from the research community on what is needed to create better tools for predicting toxicities of oligo therapeutics. Informed by results from this RFI and in alignment with NCATS’ vision of bringing more treatments for all people more quickly, we are holding this challenge competition to incentivize the development and release of publicly accessible open datasets that can be used to address bottlenecks in evaluating toxicity during drug discovery and development of oligo therapeutics.

The Challenge

To be most effective, AI-based models to predict oligo toxicity require large amounts of high-quality oligo toxicology data. To date the majority of such data is based on animal studies. Going forward, reducing the role of animal studies to predict oligo toxicity in humans would benefit from studies using in vitro human-based systems for assessing toxicology. There may also be a need to bridge the differences between predictions that are primarily based on data from animal-based studies to data collected by in vitro human-based systems. At present, however, the largest oligo toxicology datasets and predictive screening models for oligo toxicity are based on proprietary data that are not available in the public domain. These considerations highlight an overlapping need for publicly accessible high quality, open datasets that can expand use of data from in vitro human-based systems and for improving the transparency for how predictive models are developed and used in drug discovery and development.

The Oligonucleotide Toxicity (OligoTox) Open Data Challenge is a two-phase prize competition that has the goal of creating publicly accessible, high-quality datasets from in vitro human-based systems that can lead to development of improved in silico models for predicting toxicity of an oligo based on its sequence, chemical modifications, etc. The first phase (Ideation Phase) will focus on ideation by proposing data from in vitro human-based systems that can be collected and used in predicting the toxicity of a candidate oligo therapeutic. This in vitro human-based data may be supplemented by data from in vitro and in vivo animal models to improve predictive parameters and inferences. The second Data Generation Phase will incentivize the collection, generation, and contribution of this data in creating a high-quality open dataset that can be released for public access for use in training and improving predictive models of oligo toxicity. Both phases require following proper data management practices (e.g., as noted on this NIH Scientific Data Sharing webpage) as transparency and reproducibility will be critical factors for successful dissemination of these datasets for use by the broader scientific community.

In this prize challenge, NCATS seeks to incentivize the generation of high-quality and publicly accessible open datasets that make use of in vitro human systems and can be used to advance predictive models for oligo toxicity. Toxicities of interest include hepatotoxicity, kidney toxicity, thrombocytopenia, complement activation, coagulopathy, immunotoxicity, chronic neurotoxicity, and hydrocephalus. Given the availability of large datasets focused on acute neurotoxicity, specifically alterations of neuronal electrical activity, submissions focused on this topic will be considered a lower priority than other toxicities of interest.

Phase 1 Ideation
The first phase of the challenge is the Ideation Phase, and submitters are asked to describe proposed datasets that can be used to support the development and training of in silico models that can accurately predict toxicity of a candidate oligo therapeutic. These data must address one or more indicators of oligo toxicity (see list above), along with any predictor variables (e.g., oligo sequence, chemical modifications, etc.) that would be used as input features in creating computational models for predicting oligo toxicity.

Submissions should include the following:

The type(s) of oligo toxicity to be addressed
A brief summary of publicly available datasets
A discussion of how the proposed dataset will significantly add to data available in the public domain, including published scientific literature, and be useful for creating in silico predictive models using AI
A description of how the proposed dataset would be collected, generated, or provided (e.g., the intended source(s) for the dataset) and the suitability of this dataset for use in this challenge (e.g., data ownership, data use or data sharing restrictions, etc.)
A definition of the indicators and predictor variables (e.g., oligo sequence, chemical modifications, etc.) used for evaluating oligo toxicity, and evidence/reference(s) supporting how these indicators would reliably and validly predict toxicity in relation to these predictor variables
Positive and negative controls to be included, using data from the literature or other public sources
A description of the new dataset to be generated, which should include:
Model system(s) to be used, i.e., in vitro (cell culture, microphysiological system, 3-D human organoids)
Toxicology readout(s) to be assessed (e.g., cell death, changes in gene expression, etc.)
Estimated number of oligos to be evaluated, including different chemical modifications
Positive and negative control oligos and the rationale for selecting them
Methods to be used for purifying oligos, and assessing identity and purity (e.g., mass spec, HPLC, etc.)
Oligo concentrations (in vitro assays) or doses (in vivo experiments)
Number of replicates
Phase 1 submissions should be no more than 10 pages total, excluding references. The reference list should be appended separately.

In addition to a cash prize, winners of Phase 1 will be encouraged to participate in Phase 2 of the OligoTox Open Data Challenge.

Phase 2 Data Generation Phase
The second phase of this challenge is the Data Generation Phase and submitters are asked to create a high-quality, AI-ready, open dataset through the collection, generation, and contribution of data that can be used in training and improving predictive models of oligo toxicity. Only open datasets that are or will be made publicly available will be considered for a prize. Datasets based on in vitro human systems or able to extrapolate data between in vitro human systems and animal data are of particular interest. Participation in Phase 1 is not a requirement for participation in Phase 2.

For Phase 2, submitters will have 9 months to collect, build, or contribute a dataset based upon their winning submissions from Phase 1 or a dataset developed without having participated in Phase 1.

The Phase 2 submission must include four parts: a narrative document, a methodology document, a public access and dissemination plan, and a dataset.

The narrative document should include the following information, submitted as a single PDF of no more than 12 pages:

An executive summary of the dataset(s) generated, and positive/negative controls included
A summary of the main findings and conclusions
A description of how data were produced, including descriptions of relevant experimental design, data acquisition and any computational processing of raw data
A description of how indicators and predictor variables for oligo toxicity were measured, their distribution, and the distribution of predictor variables amongst tested oligos
A discussion of how the results address a gap in the publicly available data relating to one or more aspects of oligo toxicity
A discussion of how the data could be used to develop a predictive model of one or more aspects of oligo toxicity
A separate methodology document that includes information on the materials and methods used to generate the data contained in the dataset files should be submitted as a single PDF of up to 5 pages. This should include the methods used to purify and characterize oligo identity.

The dataset should consist of a data dictionary and schema documenting all metadata, and access to the raw data generated either by including a data file in Excel (or similar format) or a document with instructions for the challenge sponsors on how to access and download the raw data. This file must contain the sequences of all oligos tested, as well as the location of all chemical modifications in each oligo, data on the purity and characterization of each, and any additional metadata. Any terms for data access and data use should be defined in allowing for open and public access, such as through a creative commons license. There are no page limits for the dataset files.

Public Access and Dissemination Plan (PADP) (up to 5 pages)

Participants must submit a plan (the PADP) to maximize public access to the winning datasets from in vitro human-based systems that can lead to development of improved in silico models for predicting toxicity of an oligo based on its sequence, chemical modifications, etc.
The PADP must describe how winners will disseminate information about the solution and make the solution, as well as the knowledge necessary to access and utilize the datasets, available under non-exclusive licenses for research purposes.
The PADP must describe how the winner would allow others to utilize the solution in the event the winner is unable themselves to maximize public access to the solution. The PADP must include specific licensing schemes and scenarios when such schemes would be employed.
The PADP must describe how the winner will permit the U.S. government to allow interested parties to utilize the solution if the winner themselves fails to utilize the solution and does not permit others to utilize the solution under reasonable terms. The plan must describe any licensing schemes and scenarios that will accomplish this plan.
As described in the Rules section of this Announcement, to receive an award, Participants will be required to agree to abide by the terms of the PADP submitted by the Participants for this Challenge. Additionally, NIH intends to post or otherwise publicly display the PADP submitted by the Participants on the web or elsewhere.
Optional documentation can include documentation on the associated code used in the collection or processing of the dataset, and interactive notebooks, examples, tutorials, or documentation that demonstrate how to load and work with the dataset.

Dates:
• Phase 1 Launch: September 30, 2025
• Phase 1 Informational Webinar: January 9, 2026 (12pm-1pm ET)
• Phase 1 Submission Start/End: December 19, 2025 to February 28, 2026
• Phase 1 Judging Start/End: March 1, 2026 to April 17, 2026
• Phase 1 Winner(s) Announced: April 30, 2026
• Phase 2 Launch: May 1, 2026
• Phase 2 Informational Webinar (Tentative): May 15, 2026
• Phase 2 Submission Start/End: May 1, 2026 to December 31, 2026
• Phase 2 Judging Start/End: January 8, 2027 to February 28, 2027
• Phase 2 Winner(s) Announced: March 15, 2027

Statutory Authority to Conduct the Challenge
NCATS is conducting this Challenge under the America Creating Opportunities to Meaningfully Promote Excellence in Technology, Education, and Science (COMPETES) Reauthorization Act of 2010, as amended [15 U.S.C. § 3719]. NCATS, part of the NIH, has the mission of turning research observations into health solutions through translational science. Its statutory authority highlights the purpose of the Center to advance translational science, including by (1) coordinating and developing resources that leverage basic research in support of translational science and (2) developing partnerships and working cooperatively to foster synergy in ways that do not create duplication, redundancy, and competition with industry activities. This Challenge is consistent with and promotes NCATS’ mission by catalyzing the goal-driven development of innovative tools and technologies, e.g., oligo therapeutics, that have the potential to enhance human health.

ADDITIONAL INFORMATION
A technical assistance webinar for Phase 1 of the OligoTox Challenge will be held January 9, 2026, 12:00PM to 1:00PM Eastern Time. NIH staff will be available to answer questions related to the Challenge.

Registrants interested in the technical assistance webinar are encouraged to submit questions to ncatsoligotox@mail.nih.gov.

A technical assistance webinar for Phase 2 of the OligoTox Challenge will be held in May 2026. NIH staff will be available to answer questions related to the Challenge. Information about the webinar will be made available to registered participants.

Registrants interested in the technical assistance webinar are encouraged to submit questions by May 15, 2026 to ncatsoligotox@mail.nih.gov.

Timeline
12/19/25 09:02 AM EST: Challenge Launch (Registration opens)

12/19/25 09:03 AM EST: Phase 1 Submission Start

01/09/26 12:00 PM EST: Phase 1 Informational Webinar (Tentative)

02/28/26 11:59 PM EST: Phase 1 Submission End

03/01/26 09:00 AM EST: Phase 1 Judging Start

04/17/26 05:00 PM EDT: Phase 1 Judging End

04/30/26 05:00 PM EDT: Phase 1 Winner(s) Announced

05/01/26 09:00 AM EDT: Phase 2 Launch

05/01/26 09:00 AM EDT: Phase 2 Submission Start

06/01/26 02:00 PM ET: Phase 2 Information Webinar (Tentative)

12/31/26 11:59 PM EST: Phase 2 Submission End

01/08/27 09:00 AM EST: Phase 2 Judging Start

02/28/27 05:00 PM EST: Phase 2 Judging End

03/15/27 05:00 PM EDT: Phase 2 Winner(s) Announced

Prizes
Total cash prizes
$500,000

Amount of the Prize:
The total combined cash prize for Phases 1 and 2 is: $500,000

For Phase 1, a total prize of up to $100,000 will be awarded to the top submissions. NCATS anticipates awarding cash prizes for ten winners of up to $10,000 each.

For Phase 2, a total prize of up to $400,000 will be awarded to the top submissions. NCATS anticipates awarding cash prizes for two winners of up to $100,000 each and four runners-up with a prize of up to $50,000 each.

The number of winners may be subject to change depending on the quantity and quality of submissions. Prize funds not awarded in Phase 1 may be reallocated to Phase 2.

Award Approving Official:
The Award Approving Official will be Doctor Joni Rutter, NCATS Director.

Payment of the Prize:
Prizes awarded under this Challenge will be paid by electronic funds transfer and may be subject to federal income taxes. HHS/NIH will comply with the Internal Revenue Service withholding and reporting requirements, where applicable.

Entities participating in this Challenge are encouraged, but not required, to request and obtain a free Unique Entity ID (UEI), if they have not already done so, via SAM.gov as this will expedite prize payment. Additional information can be found at https://sam.gov/content/entity-registration.

If participating as an Individual, in the event of winning a cash prize, the Individual shall be paid the prize in full. If participating as a Team, in the event of winning a cash prize, the Team Leader shall be paid the prize in full and is solely responsible for allocating any prize amount among the members of the Team. If participating as an Entity, in the event of winning a cash prize, the prize will be paid directly to the Entity, not the Entity Point of Contact. NIH will not arbitrate, intervene, advise on, or resolve any matters among team members.

NIH reserves the right, in its sole discretion, to (a) cancel, suspend, or modify the Challenge, or any part of it, for any reason, and/or (b) not award any prizes if no submissions are deemed worthy.

Rules
Eligibility requirements
To be eligible to win a prize under this Challenge, a Participant (whether an individual, group of individuals, or entity) —

Shall have registered to participate in the Challenge under the rules promulgated by the National Institutes of Health (NIH) as published in this announcement;
Shall have complied with all the requirements set forth in this announcement;
In the case of a private entity, shall be incorporated in and maintain a primary place of business in the United States, and in the case of an individual, whether participating singly or in a group, shall be a citizen or permanent resident of the United States. However, non-U.S. citizens and non-permanent residents can participate as a member of a team that otherwise satisfies the eligibility criteria. Non-U.S. citizens and non-permanent residents are not eligible to win a monetary prize (in whole or in part). Their participation as part of a winning team, if applicable, may be recognized when the results are announced.
Shall not be a federal entity or federal employee acting within the scope of their employment;
Shall not be an employee of the Department of Health and Human Services (HHS, or any other component of HHS) acting in their personal capacity;
Who is employed by a federal agency or entity other than HHS (or any component of HHS), should consult with an agency ethics official to determine whether the federal ethics rules will limit or prohibit the acceptance of a prize under this Challenge;
Shall not be a judge of the Challenge, or any other party involved with the design, production, execution, or distribution of the Challenge or the immediate family of such a party (i.e., spouse, parent, step-parent, child, or step-child).
Shall be 18 years of age or older at the time of submission.
IP Considerations:
Through the OligoTox Challenge, NIH aims to maximize public access to the winning datasets from in vitro human-based systems that can lead to development of improved in silico models for predicting toxicity of an oligo based on its sequence, chemical modifications, etc. To balance these objectives with encouraging appropriate licensing and commercialization of the technology, the NIH is requiring participants to include a Public Access and Dissemination Plan (PADP) as part of the submission. Instructions for how to prepare this are given in the “How to Enter” section.

Public Access and Dissemination Plan (PADP):

As part of Phase 2 of the Challenge, participants were required to submit a plan (the PADP) to maximize public access to the winning datasets from in vitro human-based systems that can lead to development of improved in silico models for predicting toxicity of an oligo based on its sequence, chemical modifications, etc.
As described in the Rules section of this Announcement, to receive an award, Participants will be required to agree to abide by the terms of the PADP submitted by the Participants in Phase 2 of the Challenge. Additionally, NIH intends to post or otherwise publicly display the PADP submitted by the Participants on the web or elsewhere.

Participation Rules:

Federal grantees and recipients of cooperative agreements or other transaction (OT) awards are eligible to participate in the Challenge but may not use Federal funds from a grant award, cooperative agreement, or OT award to develop their Challenge submission or to fund efforts in support of their Challenge submission unless use of such funds is consistent with the purpose, terms, and conditions of the grant award, cooperative agreement, or OT award. Each Participant (whether participating as an individual, group of individuals, or entity) intending to use Federal grant, cooperative agreement, or OT award funds must register for and participate in the Challenge as an entity on behalf of the awardee institution, organization, or entity. If a winning Participant uses Federal grant, cooperative agreement, or OT award funds to participate in the Challenge, the prize must be treated as program income for purposes of the original grant, cooperative agreement, or OT award in accordance with applicable Uniform Administrative Requirements, Cost Principles, and Audit Requirements for Federal Awards [2 CFR § 200]. Participants using Federal grant, cooperative agreement, or OT award funds to participate and/or report prize funding as program income (for winning Participants) should coordinate with the awarding official at the federal awarding agency.
Federal contractors may not use federal funds from a contract to develop their Challenge submissions or to fund efforts in support of their Challenge submissions.
By participating in this Challenge, each Participant (whether an individual, group of individuals, or entity) agrees to assume any and all risks and waive claims against the federal government and its related entities, except in the case of willful misconduct, for any injury, death, damage, or loss of property, revenue, or profits, whether direct, indirect, or consequential, arising from participation in this Challenge, whether the injury, death, damage, or loss arises through negligence or otherwise.
Based on the subject matter of the Challenge, the type of work that it will possibly require, as well as an analysis of the likelihood of any claims for death, bodily injury, property damage, or loss potentially resulting from Challenge participation, no Participant (whether an individual, group of individuals, or entity) participating in the Challenge is required to obtain liability insurance, or demonstrate financial responsibility, or agree to indemnify the federal government against third party claims for damages arising from or related to Challenge activities in order to participate in this Challenge.
 A Participant (whether an individual, group of individuals, or entity) shall not be deemed ineligible because the Participant used federal facilities or consulted with federal employees during the Challenge if the facilities and employees are made available to all Participants participating in the Challenge on an equitable basis.
By participating in this Challenge, each Participant (whether an individual, group of individuals, or entity) warrants that they are sole author or owner of, or has the right to use, any copyrightable works that the submission comprises, that the works are wholly original with the Participant (or is an improved version of an existing work that the Participant has sufficient rights to use and improve), and that the submission does not infringe any copyright or any other rights of any third party of which the Participant is aware.
By participating in this Challenge, each Participant (whether a group of individuals or entity) grants to the NIH an irrevocable, paid-up, royalty-free nonexclusive worldwide license to reproduce, publish, post, link to, share, and display publicly the submission, including the Public Access and Dissemination Plan (PADP), submitted by the Participants, on the web or elsewhere, and agrees to grant a nonexclusive, nontransferable, irrevocable, paid-up license to practice, or have practiced for or on its behalf, the solution throughout the world. Each Participant will retain all other intellectual property rights in their submissions, as applicable. To participate in the Challenge, each Participant must warrant that there are no legal obstacles to providing the above-referenced nonexclusive licenses of the Participant’s rights to the federal government. To receive an award, Participants will not be required to transfer their intellectual property rights to NIH, but Participants must grant to the federal government the nonexclusive licenses recited herein and any license in the Participant’s PADP. As a condition for winning a cash prize in Phase 2, each Participant that has been selected as a winner will be required to agree to abide by the terms of the PADP submitted by the Participant for this Challenge.
Each Participant (whether an individual, group of individuals, or entity) agrees to follow all applicable federal, state, and local laws, regulations, and policies.
Each Participant (whether an individual, group of individuals, or entity) participating in this Challenge must comply with all terms and conditions of these rules, and participation in this Challenge constitutes each such Participant’s full and unconditional agreement to abide by these rules. Winning is contingent upon fulfilling all requirements herein.
As a condition for winning a cash prize in this Challenge, each Participant (whether an individual, group of individuals, or entity) that has been selected as a winner must complete and submit all requested winner verification and payment documents to NIH within 10 business days of formal notification. Failure to return all required verification documents by the date specified in the notification may be a basis for disqualification of a cash prize winning submission.
Judging
Phase 1 Ideation (open until 02/28/26)
Phase 2 Data Generation Phase (opens on 05/01/26)
How to enter
Phase 1 Ideation (open until 02/28/26)
Phase 2 Data Generation Phase (opens on 05/01/26)
Contact
For Further Information Contact: NCATS Challenge Prize Competitions staff at ncatsoligotox@mail.nih.gov.

Winners
Phase 1 Winners

The NIH Oligonucleotide Toxicity (OligoTox) Open Data Challenge seeks to incentivize publicly accessible high-quality datasets from human cells to predict oligonucleotide toxicity from sequence and chemical modification. Phase 1 winners have now been selected, congratulations!

Akili Bio, Inc.
Single-Cell Oligonucleotide Toxicity (SCOT) Atlas: An Extensible Platform for Profiling RNA-Induced Hepatotoxicity in Human iPSC-Derived Liver Models
Lead Innovator: Stanley Hsieh

Anote, Inc.
OligoTox Open Data Challenge – Phase 1 Ideation Submission
Lead Innovator: Natan Vidra

GOG Digital Lab
Predictive Data Frameworks for Oligonucleotide-Induced Toxicities
Lead Innovator: Gustavo Canteros

Joyce Kuo
A Scalable Multi-Task Dataset for Deep Learning Prediction of Oligonucleotide Toxicity

Lagomics Inc
Generation of a Chemistry-Resolved Interaction Dataset for AI-Based Prediction of Oligonucleotide-Induced Coagulopathy
Lead Innovator: Yasamin Jodat

LivMira Therapeutics Inc.
HuLiverTox: A Bioprinted Human Liver MPS Dataset for AI-Driven Oligonucleotide Safety Modeling
Lead Innovator: Dileep Radhamony Amma

Sentinel Gradient LLC
OligoMPS-Tox — A Multi-Organ, AI-Ready In Vitro Human Systems Dataset for Predicting Oligonucleotide Toxicity
Lead Innovator: Bhavya Johar

Uri Kartoun
Design of an AI-Ready Dataset Framework for Predicting Oligonucleotide-Induced Hepatotoxicity


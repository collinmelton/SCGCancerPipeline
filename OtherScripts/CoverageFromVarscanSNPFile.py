'''
Created on Aug 30, 2014

@author: cmelton
'''
import subprocess

def runWithUnzip(resultsFolder, outputDir, id):
#     zipfiles=map(lambda x: os.path.join(resultsFolder, id+"-", id+"-.-"+x+".varscan.snp.zip"), ["10-15","12-13","14-11","1","2","3-22","4-21","5-19","6-20","7-18","8-17","9-16","x-y-mt"])
#     for f in zipfiles: subprocess.call('unzip '+f, shell=True)
    cancerCoverage=[]
    normalCoverage=[]
    for x in ["10-15","12-13","14-11","1","2","3-22","4-21","5-19","6-20","7-18","8-17","9-16","x-y-mt"]:
        zipfile = os.path.join(resultsFolder, id+"-", id+"-.-"+x+".varscan.snp.zip")
        subprocess.call('unzip -o -d'+os.path.join(resultsFolder, id+"- ")+zipfile, shell=True)
        inputfilename = os.path.join(resultsFolder,id+"-", "mnt",id+"-disk-"+x+"-cancer-pup", id+"-.varscan.snp")
        f=open(inputfilename, "r")
        line=f.readline()
        while line[0]=="#":
            line=f.readline()
        header=line.strip().split("\t")
        normal_counts_indices=map(header.index, ["normal_reads1","normal_reads2"])
        cancer_counts_indices=map(header.index, ["tumor_reads1","tumor_reads2"])
        line=f.readline()
        while line != "":
            vals=line.strip().split("\t")
            cancerCoverage.append(sum(map(lambda x: int(vals[x]), cancer_counts_indices)))
            normalCoverage.append(sum(map(lambda x: int(vals[x]), normal_counts_indices)))
            line=f.readline()
        cancerCoverage=sorted(cancerCoverage)
        normalCoverage=sorted(normalCoverage)
        f.close()
    outputfilename=os.path.join(outputDir, "new.varscan.somatic.snp.coverage.summary")
    f=open(outputfilename, "a")
#     f.write("\t".join(["sample\taverage\tmax\tmin\tmedian"]))
    f.write("\n"+"\t".join(map(str, [id, "cancer", sum(cancerCoverage)/len(cancerCoverage), max(cancerCoverage), min(cancerCoverage), cancerCoverage[int(len(cancerCoverage)/2)]])))
    f.write("\n"+"\t".join(map(str, [id, "normal", sum(normalCoverage)/len(normalCoverage), max(normalCoverage), min(normalCoverage), normalCoverage[int(len(cancerCoverage)/2)]])))
    f.close()

def run(inputfilename, outputfilename, id):
    f=open(inputfilename, "r")
    line=f.readline()
    while line[0]=="#":
        line=f.readline()
    header=line.strip().split("\t")
    normal_counts_indices=map(header.index, ["normal_reads1","normal_reads2"])
    cancer_counts_indices=map(header.index, ["tumor_reads1","tumor_reads2"])
    cancerCoverage=[]
    normalCoverage=[]
    line=f.readline()
    while line != "":
        vals=line.strip().split("\t")
        cancerCoverage.append(sum(map(lambda x: int(vals[x]), cancer_counts_indices)))
        normalCoverage.append(sum(map(lambda x: int(vals[x]), normal_counts_indices)))
        line=f.readline()
    cancerCoverage=sorted(cancerCoverage)
    normalCoverage=sorted(normalCoverage)
    f.close()
    f=open(outputfilename, "a")
#     f.write("\t".join(["sample\taverage\tmax\tmin\tmedian"]))
    f.write("\n"+"\t".join(map(str, [id, "cancer", sum(cancerCoverage)/len(cancerCoverage), max(cancerCoverage), min(cancerCoverage), cancerCoverage[int(len(cancerCoverage)/2)]])))
    f.write("\n"+"\t".join(map(str, [id, "normal", sum(normalCoverage)/len(normalCoverage), max(normalCoverage), min(normalCoverage), normalCoverage[int(len(cancerCoverage)/2)]])))
    f.close()


varscanSNPFiles=["01eef340-598c-4205-a990-cec190ac2ca5/01eef340-598c-4205-a990-cec190ac2ca5.varscan.somatic.snp",
"028e99e9-5b9a-4954-bb6e-6d4709a3cea8/028e99e9-5b9a-4954-bb6e-6d4709a3cea8.varscan.somatic.snp",
"037c57d1-b4a5-45dc-bda4-0550461d321b/037c57d1-b4a5-45dc-bda4-0550461d321b.varscan.somatic.snp",
"045c13ef-3db7-4adf-b0a3-23338f0479f3/045c13ef-3db7-4adf-b0a3-23338f0479f3.varscan.somatic.snp",
"05fbd98c-28bd-41bc-ae60-16ac0c16c723/05fbd98c-28bd-41bc-ae60-16ac0c16c723.varscan.somatic.snp",
"084c6aec-94f7-4090-8f3f-59fa9e89721a/084c6aec-94f7-4090-8f3f-59fa9e89721a.varscan.somatic.snp",
"08da7c4c-3067-4bcf-9d7a-78566df72e69/08da7c4c-3067-4bcf-9d7a-78566df72e69.varscan.somatic.snp",
"09066fbf-6a1d-46c9-a4c5-c53340c3d1ba/09066fbf-6a1d-46c9-a4c5-c53340c3d1ba.varscan.somatic.snp",
"0956bc67-dd55-4c54-a56a-f22f2856662a/0956bc67-dd55-4c54-a56a-f22f2856662a.varscan.somatic.snp",
"0a2a3529-f645-4967-9a58-89ee20b8bb62/0a2a3529-f645-4967-9a58-89ee20b8bb62.varscan.somatic.snp",
"0ab8d063-62b4-4d47-82aa-e3351a60029d/0ab8d063-62b4-4d47-82aa-e3351a60029d.varscan.somatic.snp",
"0d66bf6c-eed0-4726-bd5b-3bf6d610b4e0/0d66bf6c-eed0-4726-bd5b-3bf6d610b4e0.varscan.somatic.snp",
"0df573ee-28f0-4244-b434-09e6ca59fbf0/0df573ee-28f0-4244-b434-09e6ca59fbf0.varscan.somatic.snp",
"0e2ee54a-51c9-4868-842d-a2a1c1cfb016/0e2ee54a-51c9-4868-842d-a2a1c1cfb016.varscan.somatic.snp",
"1174f6e4-ffbe-4e59-a000-8d861c968369/1174f6e4-ffbe-4e59-a000-8d861c968369.varscan.somatic.snp",
"11ab53c8-6366-4bc9-b60d-23e6ed2d1cae/11ab53c8-6366-4bc9-b60d-23e6ed2d1cae.varscan.somatic.snp",
"1200c8b7-a058-4c01-bbe2-04ec09a65245/1200c8b7-a058-4c01-bbe2-04ec09a65245.varscan.somatic.snp",
"13f31bcb-befd-4ed2-adc0-dc9e9ce52b75/13f31bcb-befd-4ed2-adc0-dc9e9ce52b75.varscan.somatic.snp",
"1557183c-aedd-4df6-b98a-3654667f7b69/1557183c-aedd-4df6-b98a-3654667f7b69.varscan.somatic.snp",
"15c84da3-16e5-4909-aea9-cb5894b0f8af/15c84da3-16e5-4909-aea9-cb5894b0f8af.varscan.somatic.snp",
"16e69011-c295-479f-b521-86e66fba498d/16e69011-c295-479f-b521-86e66fba498d.varscan.somatic.snp",
"180036a2-3b56-405e-a1fe-d5932517b6c7/180036a2-3b56-405e-a1fe-d5932517b6c7.varscan.somatic.snp",
"18eb4dfc-556f-4bf3-a411-4780209ed1e0/18eb4dfc-556f-4bf3-a411-4780209ed1e0.varscan.somatic.snp",
"199386c2-bb53-4fad-a1b6-59ab216a4a50/199386c2-bb53-4fad-a1b6-59ab216a4a50.varscan.somatic.snp",
"1ab7ad70-0f80-41d8-8efa-3baebc8223cc/1ab7ad70-0f80-41d8-8efa-3baebc8223cc.varscan.somatic.snp",
"1e417b93-837d-4110-9090-239eb82c3e20/1e417b93-837d-4110-9090-239eb82c3e20.varscan.somatic.snp",
"1f11de87-a5fe-49fd-80c1-e279b2bc69de/1f11de87-a5fe-49fd-80c1-e279b2bc69de.varscan.somatic.snp",
"1f744ada-6de1-446d-b62c-76a7fc8b40af/1f744ada-6de1-446d-b62c-76a7fc8b40af.varscan.somatic.snp",
"206ce0ed-36a4-47ab-ac8e-f6fca6ac5c18/206ce0ed-36a4-47ab-ac8e-f6fca6ac5c18.varscan.somatic.snp",
"20e8106b-1290-4735-abe4-7621e08e3dc8/20e8106b-1290-4735-abe4-7621e08e3dc8.varscan.somatic.snp",
"21fb46f9-4bbb-441c-af19-a687e9138344/21fb46f9-4bbb-441c-af19-a687e9138344.varscan.somatic.snp",
"220be009-8f1e-4a77-ab04-3220b5f87ea6/220be009-8f1e-4a77-ab04-3220b5f87ea6.varscan.somatic.snp",
"233b02f3-c4f0-4a67-9db5-e68d5cdaccb6/233b02f3-c4f0-4a67-9db5-e68d5cdaccb6.varscan.somatic.snp",
"267ff78b-bceb-466e-8582-c560fe227ff0/267ff78b-bceb-466e-8582-c560fe227ff0.varscan.somatic.snp",
"2779fa01-ac93-4e80-a997-3385f72172c3/2779fa01-ac93-4e80-a997-3385f72172c3.varscan.somatic.snp",
"277b02e9-ded5-4980-845d-af53690000ac/277b02e9-ded5-4980-845d-af53690000ac.varscan.somatic.snp",
"27b05b15-a44b-45ed-a6e3-e7d1ca488ea9/27b05b15-a44b-45ed-a6e3-e7d1ca488ea9.varscan.somatic.snp",
"28e9b6f9-ab6a-400b-a6a9-79095a443591/28e9b6f9-ab6a-400b-a6a9-79095a443591.varscan.somatic.snp",
"2b865f8c-539f-4cb2-8f08-41709e9511e5/2b865f8c-539f-4cb2-8f08-41709e9511e5.varscan.somatic.snp",
"2c5fa2d4-8e35-42e4-8bca-9fb3371a19c8/2c5fa2d4-8e35-42e4-8bca-9fb3371a19c8.varscan.somatic.snp",
"2c86c3ea-d926-4d39-a5ae-39ece4774287/2c86c3ea-d926-4d39-a5ae-39ece4774287.varscan.somatic.snp",
"2ed3296a-5b10-40d8-8af4-90dc031657cd/2ed3296a-5b10-40d8-8af4-90dc031657cd.varscan.somatic.snp",
"2f4a127f-101a-4192-b3e8-f9be2c8648bc/2f4a127f-101a-4192-b3e8-f9be2c8648bc.varscan.somatic.snp",
"3199cfe5-3be3-43cd-a36b-5cf2c7a9929f/3199cfe5-3be3-43cd-a36b-5cf2c7a9929f.varscan.somatic.snp",
"31c96e35-5e2f-429c-b12a-7bc5a497a300/31c96e35-5e2f-429c-b12a-7bc5a497a300.varscan.somatic.snp",
"324bcba2-f6a4-45a6-807c-215bdffcca21/324bcba2-f6a4-45a6-807c-215bdffcca21.varscan.somatic.snp",
"3258cb3b-f63b-463c-b2e4-d638149157c5/3258cb3b-f63b-463c-b2e4-d638149157c5.varscan.somatic.snp",
"328dfcb0-b5ae-43f2-bc87-dad0ab3df9e7/328dfcb0-b5ae-43f2-bc87-dad0ab3df9e7.varscan.somatic.snp",
"352768f9-3ce1-419c-beef-6515c78f5d7a/352768f9-3ce1-419c-beef-6515c78f5d7a.varscan.somatic.snp",
"35cb7841-9b09-465a-90c5-e3b8a9faad49/35cb7841-9b09-465a-90c5-e3b8a9faad49.varscan.somatic.snp",
"386b629e-fab1-4033-b088-45d6eeb4a13e/386b629e-fab1-4033-b088-45d6eeb4a13e.varscan.somatic.snp",
"38a8b734-9acc-42f9-b5b7-e51b0dfc6504/38a8b734-9acc-42f9-b5b7-e51b0dfc6504.varscan.somatic.snp",
"3a4a6aff-4720-41fc-a114-85f38610ebb1/3a4a6aff-4720-41fc-a114-85f38610ebb1.varscan.somatic.snp",
"3c975b68-4a6e-4c97-a8d1-b9a427d4957c/3c975b68-4a6e-4c97-a8d1-b9a427d4957c.varscan.somatic.snp",
"3cbca837-f5a7-4a87-8f02-c59eac232d5a/3cbca837-f5a7-4a87-8f02-c59eac232d5a.varscan.somatic.snp",
"3d1f4059-2220-45b4-a4d2-b14f76cec96a/3d1f4059-2220-45b4-a4d2-b14f76cec96a.varscan.somatic.snp",
"3e976301-4042-4ac0-8f01-ef2f1a8161a6/3e976301-4042-4ac0-8f01-ef2f1a8161a6.varscan.somatic.snp",
"3edf933a-be1a-46e8-bfb0-acdff30e64a0/3edf933a-be1a-46e8-bfb0-acdff30e64a0.varscan.somatic.snp",
"3f93d1cf-a307-4e29-9e6f-239e9ffaa179/3f93d1cf-a307-4e29-9e6f-239e9ffaa179.varscan.somatic.snp",
"3f960d3b-a58c-43d0-a8a4-f3555b399c9d/3f960d3b-a58c-43d0-a8a4-f3555b399c9d.varscan.somatic.snp",
"4075717b-2240-4842-a00a-0c4dbcb910c4/4075717b-2240-4842-a00a-0c4dbcb910c4.varscan.somatic.snp",
"43abe847-4ba7-466e-8283-5d7b80b999a7/43abe847-4ba7-466e-8283-5d7b80b999a7.varscan.somatic.snp",
"441919eb-8c9f-478c-bda3-8de4f8295e8a/441919eb-8c9f-478c-bda3-8de4f8295e8a.varscan.somatic.snp",
"446ce2a3-d328-443c-a419-3344baad0e16/446ce2a3-d328-443c-a419-3344baad0e16.varscan.somatic.snp",
"44dec838-b653-42a7-a58b-a1fd232cd68c/44dec838-b653-42a7-a58b-a1fd232cd68c.varscan.somatic.snp",
"46592b7b-6968-42a6-83af-0917c9f4a9a5/46592b7b-6968-42a6-83af-0917c9f4a9a5.varscan.somatic.snp",
"468649e7-1525-4b8a-8ab5-ddb27db5f022/468649e7-1525-4b8a-8ab5-ddb27db5f022.varscan.somatic.snp",
"4b2ec8aa-8460-440e-b75b-6b92ae7a3ffc/4b2ec8aa-8460-440e-b75b-6b92ae7a3ffc.varscan.somatic.snp",
"4bf50455-9ab7-4521-a791-089f66d3b877/4bf50455-9ab7-4521-a791-089f66d3b877.varscan.somatic.snp",
"4bfbce2b-9d0b-4e8a-950f-fd8e0ba3e05a/4bfbce2b-9d0b-4e8a-950f-fd8e0ba3e05a.varscan.somatic.snp",
"4c31127e-d095-4978-9dac-35153c27f6ed/4c31127e-d095-4978-9dac-35153c27f6ed.varscan.somatic.snp",
"4c325ee1-9b8e-45d4-a0a7-e47ec1eac402/4c325ee1-9b8e-45d4-a0a7-e47ec1eac402.varscan.somatic.snp",
"4c42dc4e-66b7-40bf-9fbe-b92543248198/4c42dc4e-66b7-40bf-9fbe-b92543248198.varscan.somatic.snp",
"4cffea0b-90a7-4c86-a73f-bb8feca3ada7/4cffea0b-90a7-4c86-a73f-bb8feca3ada7.varscan.somatic.snp",
"4edff57f-4b0e-4770-beac-590da7d7232c/4edff57f-4b0e-4770-beac-590da7d7232c.varscan.somatic.snp",
"4fca10bd-988a-47f9-8ed5-037c0d59d70c/4fca10bd-988a-47f9-8ed5-037c0d59d70c.varscan.somatic.snp",
"521ea765-1bd1-423d-a75d-091243df37a9/521ea765-1bd1-423d-a75d-091243df37a9.varscan.somatic.snp",
"5338d435-68fb-4f0d-a3e6-c843f703f75f/5338d435-68fb-4f0d-a3e6-c843f703f75f.varscan.somatic.snp",
"534f85f2-b9eb-45be-87d5-c3dffd13d541/534f85f2-b9eb-45be-87d5-c3dffd13d541.varscan.somatic.snp",
"54254f5a-50e8-4150-a9b7-56a0470d2a56/54254f5a-50e8-4150-a9b7-56a0470d2a56.varscan.somatic.snp",
"5580b21a-2cdb-4777-ad79-6e06654144f5/5580b21a-2cdb-4777-ad79-6e06654144f5.varscan.somatic.snp",
"5613f5f6-086d-470b-8f77-6dcb7f8625b7/5613f5f6-086d-470b-8f77-6dcb7f8625b7.varscan.somatic.snp",
"566792ae-f853-4a47-856d-f02cdcfcb18a/566792ae-f853-4a47-856d-f02cdcfcb18a.varscan.somatic.snp",
"585e6487-b0a3-4828-8a06-46bee01dff74/585e6487-b0a3-4828-8a06-46bee01dff74.varscan.somatic.snp",
"5a57dc25-d252-4b22-b192-4de630d7002f/5a57dc25-d252-4b22-b192-4de630d7002f.varscan.somatic.snp",
"61c655ec-52b5-453f-a6cc-b2aba445b027/61c655ec-52b5-453f-a6cc-b2aba445b027.varscan.somatic.snp",
"63bd2175-4b7c-44c9-aef3-9efc8f79837b/63bd2175-4b7c-44c9-aef3-9efc8f79837b.varscan.somatic.snp",
"659668e8-f0d9-4ff2-bbc8-9246f2ef49ab/659668e8-f0d9-4ff2-bbc8-9246f2ef49ab.varscan.somatic.snp",
"673493f6-975c-49e8-934c-001e9a0fff90/673493f6-975c-49e8-934c-001e9a0fff90.varscan.somatic.snp",
"6956473b-92e5-4069-a26b-b48e280e76f2/6956473b-92e5-4069-a26b-b48e280e76f2.varscan.somatic.snp",
"69d56f2d-6924-409b-9d1e-c8d69b400270/69d56f2d-6924-409b-9d1e-c8d69b400270.varscan.somatic.snp",
"6c28f086-6a25-40b6-93eb-bba0014acda6/6c28f086-6a25-40b6-93eb-bba0014acda6.varscan.somatic.snp",
"6c5154d2-af36-492f-b520-d925528824e4/6c5154d2-af36-492f-b520-d925528824e4.varscan.somatic.snp",
"6d72de06-232a-4983-a06c-eba6d82cb3f1/6d72de06-232a-4983-a06c-eba6d82cb3f1.varscan.somatic.snp",
"6dfd47d2-831a-4386-9051-f78199a16bb5/6dfd47d2-831a-4386-9051-f78199a16bb5.varscan.somatic.snp",
"6fd72426-f6c8-47ca-a500-d5d3600b9b15/6fd72426-f6c8-47ca-a500-d5d3600b9b15.varscan.somatic.snp",
"713190ed-c6c1-4695-814b-85ca9b95e6a0/713190ed-c6c1-4695-814b-85ca9b95e6a0.varscan.somatic.snp",
"71fed998-8fac-4e6b-b38a-22c83e05e958/71fed998-8fac-4e6b-b38a-22c83e05e958.varscan.somatic.snp",
"733a2a0f-b37a-4b81-b49e-3c0f30d1eb37/733a2a0f-b37a-4b81-b49e-3c0f30d1eb37.varscan.somatic.snp",
"737b35e1-d668-4fce-9b6e-76946c7952b6/737b35e1-d668-4fce-9b6e-76946c7952b6.varscan.somatic.snp",
"762dea8a-5b41-4058-979a-b7876ed13d7e/762dea8a-5b41-4058-979a-b7876ed13d7e.varscan.somatic.snp",
"781acafa-a58b-4180-bbca-c5b675d51b17/781acafa-a58b-4180-bbca-c5b675d51b17.varscan.somatic.snp",
"784de7ac-8424-42eb-83d4-a1bebaa42b97/784de7ac-8424-42eb-83d4-a1bebaa42b97.varscan.somatic.snp",
"7861a5f2-8910-4bcb-9c34-f79c5acd6e21/7861a5f2-8910-4bcb-9c34-f79c5acd6e21.varscan.somatic.snp",
"786e8dbe-442e-4551-87b3-b4c333b04dd4/786e8dbe-442e-4551-87b3-b4c333b04dd4.varscan.somatic.snp",
"7b0622ab-63ea-483f-ae40-d3ea587bdbba/7b0622ab-63ea-483f-ae40-d3ea587bdbba.varscan.somatic.snp",
"7b7d2d2b-abe6-4f32-add0-3715096995d1/7b7d2d2b-abe6-4f32-add0-3715096995d1.varscan.somatic.snp",
"7b9947d5-e86b-4006-8890-ea5209584a88/7b9947d5-e86b-4006-8890-ea5209584a88.varscan.somatic.snp",
"7d3e77cc-9603-4722-8e39-1912b678871b/7d3e77cc-9603-4722-8e39-1912b678871b.varscan.somatic.snp",
"7dce4215-01f6-43e9-9378-9e22c26c23a2/7dce4215-01f6-43e9-9378-9e22c26c23a2.varscan.somatic.snp",
"7de19081-d5fd-468c-ad0d-f6e3e8b2ad70/7de19081-d5fd-468c-ad0d-f6e3e8b2ad70.varscan.somatic.snp",
"80c4fc5b-dc7d-4308-8c81-6e83625963b5/80c4fc5b-dc7d-4308-8c81-6e83625963b5.varscan.somatic.snp",
"832316ee-08ab-412b-9e4d-8a57ff84bc11/832316ee-08ab-412b-9e4d-8a57ff84bc11.varscan.somatic.snp",
"8994c8c4-9ea0-46cb-b4a5-0055d7da5bfa/8994c8c4-9ea0-46cb-b4a5-0055d7da5bfa.varscan.somatic.snp",
"89f91e36-fed8-4da3-8a35-c761d6f65285/89f91e36-fed8-4da3-8a35-c761d6f65285.varscan.somatic.snp",
"8aaa4e25-5c12-4ace-96dc-91aaa0c4457c/8aaa4e25-5c12-4ace-96dc-91aaa0c4457c.varscan.somatic.snp",
"8bd68696-5eeb-46d7-9dc6-c39c849029a4/8bd68696-5eeb-46d7-9dc6-c39c849029a4.varscan.somatic.snp",
"8cf54607-01ce-42b0-9bd9-8627edd9f3b7/8cf54607-01ce-42b0-9bd9-8627edd9f3b7.varscan.somatic.snp",
"900b5f21-e8fb-4a3f-b8e3-ef098716e508/900b5f21-e8fb-4a3f-b8e3-ef098716e508.varscan.somatic.snp",
"90a8d2cb-1f31-45f9-8575-eed3cb9a7798/90a8d2cb-1f31-45f9-8575-eed3cb9a7798.varscan.somatic.snp",
"914d8613-3c25-4a10-be50-28b42f4d3a5d/914d8613-3c25-4a10-be50-28b42f4d3a5d.varscan.somatic.snp",
"9293e197-e38a-4e19-a7d0-1b45d1ad48bd/9293e197-e38a-4e19-a7d0-1b45d1ad48bd.varscan.somatic.snp",
"93ed7a2b-b0cb-4a84-871f-5c34a0b6a640/93ed7a2b-b0cb-4a84-871f-5c34a0b6a640.varscan.somatic.snp",
"9435447e-d65f-408b-863b-6576b1d652dd/9435447e-d65f-408b-863b-6576b1d652dd.varscan.somatic.snp",
"96415324-6859-41eb-8d07-be4e6f7d9781/96415324-6859-41eb-8d07-be4e6f7d9781.varscan.somatic.snp",
"9853c6bb-a42c-4698-abe5-f3c897ebc8f6/9853c6bb-a42c-4698-abe5-f3c897ebc8f6.varscan.somatic.snp",
"993103b1-e5a1-4c33-8629-be53ebc41d64/993103b1-e5a1-4c33-8629-be53ebc41d64.varscan.somatic.snp",
"9938ce5c-e74e-446e-a932-f096f85cc3b1/9938ce5c-e74e-446e-a932-f096f85cc3b1.varscan.somatic.snp",
"99e32c59-aa73-43fa-88c9-399ddadb2c72/99e32c59-aa73-43fa-88c9-399ddadb2c72.varscan.somatic.snp",
"9aa8c2b4-db32-4899-bf0c-d578cb175e90/9aa8c2b4-db32-4899-bf0c-d578cb175e90.varscan.somatic.snp",
"9ae5101b-6031-4ab1-bfca-198c41af3184/9ae5101b-6031-4ab1-bfca-198c41af3184.varscan.somatic.snp",
"9b84d013-713b-4571-888c-12863fdd4a52/9b84d013-713b-4571-888c-12863fdd4a52.varscan.somatic.snp",
"9d95a65b-e41d-4f93-92d4-99dce29ff40d/9d95a65b-e41d-4f93-92d4-99dce29ff40d.varscan.somatic.snp",
"9da334ed-314c-4445-ba29-ffd9bb42403c/9da334ed-314c-4445-ba29-ffd9bb42403c.varscan.somatic.snp",
"9dc7812b-c7a2-4de4-bf6d-4c7261384a62/9dc7812b-c7a2-4de4-bf6d-4c7261384a62.varscan.somatic.snp",
"9e178d80-fe90-4998-9dae-726c4beceef7/9e178d80-fe90-4998-9dae-726c4beceef7.varscan.somatic.snp",
"9f5c0a43-d222-484a-9249-d288387dd7a1/9f5c0a43-d222-484a-9249-d288387dd7a1.varscan.somatic.snp",
"9f632fc3-2241-469a-a14a-cd537350155f/9f632fc3-2241-469a-a14a-cd537350155f.varscan.somatic.snp",
"9f81c602-8afa-4588-b0b6-6e5a1a128d5a/9f81c602-8afa-4588-b0b6-6e5a1a128d5a.varscan.somatic.snp",
"9f89510c-ed07-471f-b35e-7c87c237b9fe/9f89510c-ed07-471f-b35e-7c87c237b9fe.varscan.somatic.snp",
"a02932f3-5606-42a5-aee4-e5b9f5aa1aaa/a02932f3-5606-42a5-aee4-e5b9f5aa1aaa.varscan.somatic.snp",
"a081da6f-51f3-4555-a6b8-736f8c6cda39/a081da6f-51f3-4555-a6b8-736f8c6cda39.varscan.somatic.snp",
"a1e65587-24c1-4b41-92a7-4e1f15fffd78/a1e65587-24c1-4b41-92a7-4e1f15fffd78.varscan.somatic.snp",
"a2ac9937-f351-4d78-9261-264bf6c21e0c/a2ac9937-f351-4d78-9261-264bf6c21e0c.varscan.somatic.snp",
"a3de401d-91fe-49a2-bb07-81c1a06506e6/a3de401d-91fe-49a2-bb07-81c1a06506e6.varscan.somatic.snp",
"a3e1ac67-a1f2-44fb-8343-a7e8239fc24a/a3e1ac67-a1f2-44fb-8343-a7e8239fc24a.varscan.somatic.snp",
"a515cf2d-e918-4958-9bf6-e611b425a97e/a515cf2d-e918-4958-9bf6-e611b425a97e.varscan.somatic.snp",
"a587e62e-430e-4f1b-82b7-6fbe856fdaf1/a587e62e-430e-4f1b-82b7-6fbe856fdaf1.varscan.somatic.snp",
"a8015490-9740-45c9-8bd2-eb6d1beefc2e/a8015490-9740-45c9-8bd2-eb6d1beefc2e.varscan.somatic.snp",
"a824b3bd-34d5-4cc1-a92f-f9d6ac0f1814/a824b3bd-34d5-4cc1-a92f-f9d6ac0f1814.varscan.somatic.snp",
"a8614c49-c7f1-4de9-bd2f-4764b094d604/a8614c49-c7f1-4de9-bd2f-4764b094d604.varscan.somatic.snp",
"a8d6694c-a213-4544-ac0b-63bce16d8f4e/a8d6694c-a213-4544-ac0b-63bce16d8f4e.varscan.somatic.snp",
"a9b7d7fe-be31-4f71-afee-c1bfdf511888/a9b7d7fe-be31-4f71-afee-c1bfdf511888.varscan.somatic.snp",
"a9d3a7b0-faf8-4e56-8600-d9e6882a4f23/a9d3a7b0-faf8-4e56-8600-d9e6882a4f23.varscan.somatic.snp",
"acae4e48-57de-40c5-9ae5-4386144ebaea/acae4e48-57de-40c5-9ae5-4386144ebaea.varscan.somatic.snp",
"af37d8f0-fea0-47f6-be0a-2fb55a5da5c4/af37d8f0-fea0-47f6-be0a-2fb55a5da5c4.varscan.somatic.snp",
"b1c4579c-72fa-4974-830c-b6a81ace45f4/b1c4579c-72fa-4974-830c-b6a81ace45f4.varscan.somatic.snp",
"b3631718-9e0a-454c-bee1-8f36ebc509d8/b3631718-9e0a-454c-bee1-8f36ebc509d8.varscan.somatic.snp",
"b37d944c-c65c-46d6-bb8b-3eb37cb85b68/b37d944c-c65c-46d6-bb8b-3eb37cb85b68.varscan.somatic.snp",
"b43e41af-1d82-4b5f-b8f1-add0510e6b86/b43e41af-1d82-4b5f-b8f1-add0510e6b86.varscan.somatic.snp",
"b58ad350-5140-4fa8-bc2c-24bca8395f3a/b58ad350-5140-4fa8-bc2c-24bca8395f3a.varscan.somatic.snp",
"b5fba77b-1f50-4e71-95b2-566afba4bdd7/b5fba77b-1f50-4e71-95b2-566afba4bdd7.varscan.somatic.snp",
"b60f22ac-a659-4f33-b01d-820e86a9a5c9/b60f22ac-a659-4f33-b01d-820e86a9a5c9.varscan.somatic.snp",
"b913d254-8307-4b8a-8313-d978e32bb38f/b913d254-8307-4b8a-8313-d978e32bb38f.varscan.somatic.snp",
"ba9d4568-71a8-41fc-a2be-0bf784e36d6a/ba9d4568-71a8-41fc-a2be-0bf784e36d6a.varscan.somatic.snp",
"bbdaa931-e922-49be-bfbf-fa0c2ae27d7a/bbdaa931-e922-49be-bfbf-fa0c2ae27d7a.varscan.somatic.snp",
"bd4284a7-260d-48b3-9b48-a1c578c30d27/bd4284a7-260d-48b3-9b48-a1c578c30d27.varscan.somatic.snp",
"bd75d8ee-916b-4abc-bf33-bce6a4217076/bd75d8ee-916b-4abc-bf33-bce6a4217076.varscan.somatic.snp",
"be3a7ef3-34ed-40e1-9d9c-187940596b26/be3a7ef3-34ed-40e1-9d9c-187940596b26.varscan.somatic.snp",
"bf15f7ad-9d92-473b-91d1-f24aa373ab97/bf15f7ad-9d92-473b-91d1-f24aa373ab97.varscan.somatic.snp",
"bf6b931e-cb2e-4cc0-b603-e9fc71fdc509/bf6b931e-cb2e-4cc0-b603-e9fc71fdc509.varscan.somatic.snp",
"bff84539-7862-45b2-b5fc-e77291fcca8b/bff84539-7862-45b2-b5fc-e77291fcca8b.varscan.somatic.snp",
"c0892598-1f7b-4f23-9cd8-731f797753d5/c0892598-1f7b-4f23-9cd8-731f797753d5.varscan.somatic.snp",
"c184c3ca-7ad3-4202-b108-cb9fd5f5d947/c184c3ca-7ad3-4202-b108-cb9fd5f5d947.varscan.somatic.snp",
"c1f50a22-38df-41cc-a1f4-f7985504a7ac/c1f50a22-38df-41cc-a1f4-f7985504a7ac.varscan.somatic.snp",
"c350b53b-3c27-4b9e-8b57-d7b6a1205a1f/c350b53b-3c27-4b9e-8b57-d7b6a1205a1f.varscan.somatic.snp",
"c364e81c-eb1e-4870-ab37-9c661f5f2e3d/c364e81c-eb1e-4870-ab37-9c661f5f2e3d.varscan.somatic.snp",
"c931f3bd-74c5-4ebd-bc0f-c7c6becd25ab/c931f3bd-74c5-4ebd-bc0f-c7c6becd25ab.varscan.somatic.snp",
"cc353686-e2b6-46a1-8095-48e2bc44b11f/cc353686-e2b6-46a1-8095-48e2bc44b11f.varscan.somatic.snp",
"cc4bd56a-25c5-4c48-b583-ac3aeb778ca6/cc4bd56a-25c5-4c48-b583-ac3aeb778ca6.varscan.somatic.snp",
"ccd4a24b-d8cc-4686-9dee-c98b0c5a8d21/ccd4a24b-d8cc-4686-9dee-c98b0c5a8d21.varscan.somatic.snp",
"cd032ddb-55f2-4c77-8dcf-e4e630f7de6f/cd032ddb-55f2-4c77-8dcf-e4e630f7de6f.varscan.somatic.snp",
"cd0dc947-b708-464e-beb7-954c9d4583e7/cd0dc947-b708-464e-beb7-954c9d4583e7.varscan.somatic.snp",
"ce8612ab-3149-4a6a-b424-29c0c21c9b8b/ce8612ab-3149-4a6a-b424-29c0c21c9b8b.varscan.somatic.snp",
"cfc999aa-802e-4a0d-8b14-be2b19fa2a17/cfc999aa-802e-4a0d-8b14-be2b19fa2a17.varscan.somatic.snp",
"d17668c4-5d66-4e16-8b70-b54616581aae/d17668c4-5d66-4e16-8b70-b54616581aae.varscan.somatic.snp",
"d189a6c6-210f-4957-b500-457d5ce3f867/d189a6c6-210f-4957-b500-457d5ce3f867.varscan.somatic.snp",
"d1eeaca3-52c0-4a0c-916a-fd64cec911d9/d1eeaca3-52c0-4a0c-916a-fd64cec911d9.varscan.somatic.snp",
"d2824e6d-3784-45c2-9b0f-52b17356b5da/d2824e6d-3784-45c2-9b0f-52b17356b5da.varscan.somatic.snp",
"d2aa7d70-a25f-4a9e-9d52-37e372461824/d2aa7d70-a25f-4a9e-9d52-37e372461824.varscan.somatic.snp",
"d2ca9a4d-c14c-40c0-871c-88794a02e4aa/d2ca9a4d-c14c-40c0-871c-88794a02e4aa.varscan.somatic.snp",
"d2f2560c-ec80-4fea-9474-c47a2e85ea95/d2f2560c-ec80-4fea-9474-c47a2e85ea95.varscan.somatic.snp",
"d5326429-9805-47f9-97b0-fbda658e3f01/d5326429-9805-47f9-97b0-fbda658e3f01.varscan.somatic.snp",
"d5eedfb5-4fa3-4a1e-ad0d-ea550b22a5ba/d5eedfb5-4fa3-4a1e-ad0d-ea550b22a5ba.varscan.somatic.snp",
"d67cd793-2931-429a-9084-2f3c4c8be7ad/d67cd793-2931-429a-9084-2f3c4c8be7ad.varscan.somatic.snp",
"d7641ca1-a062-4dd2-9414-447af7311a84/d7641ca1-a062-4dd2-9414-447af7311a84.varscan.somatic.snp",
"d8d7a6c2-6427-4f47-968c-6c3affba4617/d8d7a6c2-6427-4f47-968c-6c3affba4617.varscan.somatic.snp",
"da61857a-b00e-4a19-ab02-8fbb2ee56bea/da61857a-b00e-4a19-ab02-8fbb2ee56bea.varscan.somatic.snp",
"db058eb4-9c5e-4226-92d1-de24aea83630/db058eb4-9c5e-4226-92d1-de24aea83630.varscan.somatic.snp",
"dbd5b0de-94c9-45dd-afb3-6820a7ecaca2/dbd5b0de-94c9-45dd-afb3-6820a7ecaca2.varscan.somatic.snp",
"dcff8659-ecca-4798-80a3-dd6e7642ba2b/dcff8659-ecca-4798-80a3-dd6e7642ba2b.varscan.somatic.snp",
"dee43082-4aa7-4650-981c-a779e9fe7b55/dee43082-4aa7-4650-981c-a779e9fe7b55.varscan.somatic.snp",
"e16ca88f-488b-40f0-9169-e5a62482a2ff/e16ca88f-488b-40f0-9169-e5a62482a2ff.varscan.somatic.snp",
"e1e8dc66-467c-48a5-be56-3db4f665b9f1/e1e8dc66-467c-48a5-be56-3db4f665b9f1.varscan.somatic.snp",
"e4fc0909-f284-4471-866d-d8967b6adcbc/e4fc0909-f284-4471-866d-d8967b6adcbc.varscan.somatic.snp",
"e58e8850-154f-4695-bee0-005c76410327/e58e8850-154f-4695-bee0-005c76410327.varscan.somatic.snp",
"e737f650-b72d-44e7-b750-558a56716803/e737f650-b72d-44e7-b750-558a56716803.varscan.somatic.snp",
"e9f4f373-37a5-48ad-a1a0-b0d47820111a/e9f4f373-37a5-48ad-a1a0-b0d47820111a.varscan.somatic.snp",
"ea3b5da2-6a12-400c-bf0f-e442f5ec1132/ea3b5da2-6a12-400c-bf0f-e442f5ec1132.varscan.somatic.snp",
"eb2dbb4f-66b6-4525-8323-431970f7a64e/eb2dbb4f-66b6-4525-8323-431970f7a64e.varscan.somatic.snp",
"ec7a9d74-8b3f-4747-8926-1739a016ab2b/ec7a9d74-8b3f-4747-8926-1739a016ab2b.varscan.somatic.snp",
"ee0a4a13-613e-4c5d-96c3-8083a013702d/ee0a4a13-613e-4c5d-96c3-8083a013702d.varscan.somatic.snp",
"ef6f6553-2575-457e-bc07-401215e54759/ef6f6553-2575-457e-bc07-401215e54759.varscan.somatic.snp",
"f1985735-a188-4567-a88a-530f9b80291b/f1985735-a188-4567-a88a-530f9b80291b.varscan.somatic.snp",
"f2801b21-5444-4cc3-a642-c60c6d82cd3d/f2801b21-5444-4cc3-a642-c60c6d82cd3d.varscan.somatic.snp",
"f34aa3b6-e966-49c6-bc55-130545772c53/f34aa3b6-e966-49c6-bc55-130545772c53.varscan.somatic.snp",
"f45f4a30-6b4c-4f15-9140-959d6a25a45f/f45f4a30-6b4c-4f15-9140-959d6a25a45f.varscan.somatic.snp",
"f5759059-c0e3-4f1a-af96-5c7197d3c33c/f5759059-c0e3-4f1a-af96-5c7197d3c33c.varscan.somatic.snp",
"f5bc5d97-e054-4e53-992d-71b896bd97d5/f5bc5d97-e054-4e53-992d-71b896bd97d5.varscan.somatic.snp",
"f6e1ec78-5ad9-4879-9b7e-262d17b166ad/f6e1ec78-5ad9-4879-9b7e-262d17b166ad.varscan.somatic.snp",
"f77b6930-a1ce-446e-a15d-c018cbbecfee/f77b6930-a1ce-446e-a15d-c018cbbecfee.varscan.somatic.snp",
"f77cf4d3-9829-4fd3-95a1-b29aed9e94b8/f77cf4d3-9829-4fd3-95a1-b29aed9e94b8.varscan.somatic.snp",
"f8ac5fd2-11f2-4b3a-a3a6-ff1a6d8e1cae/f8ac5fd2-11f2-4b3a-a3a6-ff1a6d8e1cae.varscan.somatic.snp",
"fb79c491-7b01-42ae-8369-8364e442e31b/fb79c491-7b01-42ae-8369-8364e442e31b.varscan.somatic.snp",
"fba80122-d8b2-4d8d-a032-9767e8160f9f/fba80122-d8b2-4d8d-a032-9767e8160f9f.varscan.somatic.snp",
"fcecafa8-8b8b-4ae6-81bc-73acf05b211a/fcecafa8-8b8b-4ae6-81bc-73acf05b211a.varscan.somatic.snp",
"fd9ee494-65fe-4de4-adff-7952a059b17f/fd9ee494-65fe-4de4-adff-7952a059b17f.varscan.somatic.snp",
"ff530f28-0ec0-4494-bb54-44bb055bae1c/ff530f28-0ec0-4494-bb54-44bb055bae1c.varscan.somatic.snp",
"ff9def3d-17e5-4ef6-b74e-933f11ed6f00/ff9def3d-17e5-4ef6-b74e-933f11ed6f00.varscan.somatic.snp",
"ffaa98a0-2b69-46dc-aee5-c5c3f2abbc38/ffaa98a0-2b69-46dc-aee5-c5c3f2abbc38.varscan.somatic.snp"]

import os

# WGSDir="/srv/gsfs0/projects/snyder/collinmelton/TCGA_Variants/WGSResults"
WGSDir="/srv/gsfs0/projects/gbsc/Clinical/TCGA_Variants/gbsc-gcp-lab-snyder-users-cmelton/results"
outputDir="/srv/gsfs0/projects/snyder/collinmelton/TCGA_Variants/WGSCoverage/"
#complete "id797a71ee-372d-41e5-aeee-5ab3c4661110","id799e9106-e62e-4681-8861-137f474ee8b8","id0078b0c4-68a9-483b-9aab-61156d263213","id79ecc142-8f6e-4260-a7ce-c8e2a6212c35","id7ad0617e-bc4a-4869-8419-744f04d9b7a3","id00ad0ffe-2105-4829-a495-1c2aceb5bb31","id7b982d5e-3a7d-40ac-bd25-6044c62879b6","id00bca18c-b3d4-45a3-8f19-034cc40449a4","id7c93642d-7e05-40f8-b1ef-a014edcfba42","id010aac75-3bfe-4bf2-b866-af0f2d92f125","id7ca06dbe-f073-4aee-86ff-a545fad5ea93","id01a92062-967a-4900-8dc7-a5ecd3b3f8e2","id7ca97692-0bf5-4bbd-81ce-10a051d04bd5","id01c4ca12-ac08-41dd-9e4f-bfb0971688a2","id7d905f8f-3967-4ea8-96c8-17b1b03fbec3",
pids=["id7dcc809b-e33a-4453-b92a-c00786f48cb0","id7dd3bda9-3a0a-46fa-b052-a6d24b847804","id7fee9d6f-c0ea-43c2-89c2-66d965a077a8","id81b7cbc1-c037-426f-95b0-d729a30697da","id827453a9-7388-4336-9ac5-4defcb904607","id830fd359-abd4-4c99-89e4-2d6bb9ce5637","id83b3060c-2449-4581-88ef-817d126e4525","id84569c1a-6baa-4380-ae3c-707df1be4618","id84ae73e4-3b39-4444-9bfb-3fc6057eab32","id8785012f-f73e-4d68-87cf-1d804af32782","id87ac0ad4-8c9d-409e-9a86-1d201d01769d","id88db1340-e4bf-451a-87c0-6e9168296f5e","id8a4da4f0-30f2-497e-8a7c-988fd6b813cb","id02ec7a3a-2812-4a65-93cd-07bf655cc91b","id8a6d2ce3-cc57-451b-9b07-8263782aa23f","id0398eae1-7216-4595-80a5-6b117d96e070","id8afa3140-dcbd-44c0-ab13-4ada6a0444d4","id03a9dd9d-62ae-4acd-9272-389274858f3d","id8b5746f9-dbee-40bd-9141-1081960cf286","id03f3dd32-e1d1-485b-a968-3f55798f4d46","id8c238d30-df8e-4e6b-98fc-21696269a294","id0487fc41-386c-4d76-9084-daf959bf5e98","id8cb0144b-be6b-40a1-86a2-708f96d9b615","id054f472f-98cb-4559-b2e2-b5f800fc8eef","id8d3e3d09-bc61-4d5e-b8d6-5aee6cb488ae","id05506f4c-e701-4a9d-ae06-97f066aade43","id8d54a2ed-03f1-4a23-bdd2-f3395f5d3716","id072507cf-1249-4b28-86d4-3ab97822a74c","id8d9e4917-334b-4c76-aee1-1e22be772db0","id0798dbe2-1914-427c-a2fe-2a865d0d6eda","id8da3103e-3e6c-4176-a583-d5fe5e60601e","id0809ba8b-4ab6-4f43-934c-c1ccbc014a7e","id8e6020a1-68df-4d57-ad4e-33acc0861fb3","id085ffcdc-8600-42e2-852a-506465f716a2","id8f3acbbd-b2ca-4e66-8fe5-dbb49277f35b","id09587f1c-5c99-4102-bc49-84d50fa8d0ce","id8fc1f1be-d2d5-4b3a-9973-f4d964018beb","id096b4f32-10c1-4737-a0dd-cae04c54ee33","id8fccff71-3710-4cb9-b8af-61707aebe5d8","id098acb76-0bf5-44e5-bcae-f919cf5fa5e5","id904a8757-e0c5-41ef-b583-c8f170caaac1","id0a48873a-092a-4b25-822b-b0b3c18c08a8","id90f44998-3c05-44b2-9501-bbd0f43c5147","id0aecac64-5982-4d76-8f31-958f6a00951d","id90f4b65c-cfd4-4066-a5b9-842885c172e2","id0bc5744c-5fa3-45bb-87d0-70a02068b392","id9205c164-7975-421a-90c3-edfe8def595c","id0c1a2e7d-e7e4-481e-a012-ef214c444497","id925ae8ff-3200-4564-be76-94db0984b027","id0c50a2c2-1d4c-45ca-9127-f96d4db18daf","id932bcdd4-5145-49e8-9a95-8a546396f572","id0dca98b0-f43e-45b6-9a02-00092c78678c","id934a58af-0ffc-487e-9e95-fa7643b7d5cb","id0eda7373-8335-455f-b90a-8c3540f28088","id93ed7a2b-b0cb-4a84-871f-5c34a0b6a640","id0f3c4f8f-0766-46b0-81b4-18ca4f6ea20a","id941f75a1-fea4-4539-ba69-60bb11608f6d","id0f530b3e-5b6d-4892-9ebd-9138d76fdca7","id944a102f-3475-45aa-9e30-4b92d1449f00","id12581634-eebb-4841-8498-71dc9d78c546","id94817e78-78f5-4556-9394-610df65605e2","id12824b2b-9c0a-4dcf-8941-34593a4e93da","id96312510-c126-485d-8109-ed81844a1dc3","id13f5814c-1f99-4ffa-84a4-3bbd8979faae","id96457ebf-ae88-41f4-97e5-cdc3163418aa","id1549dc64-3dab-43fc-96e9-b07d520957e1","id97686ab2-3b7e-4b21-9bf3-9a5a01709e04","id15556b28-c6bd-455a-87b6-4b7b3e33d0e4","id97831d28-ab41-4c18-bfc2-c4c6bc757d13","id174850b4-5ec2-462b-a890-89bd1716b3c2","id9853c6bb-a42c-4698-abe5-f3c897ebc8f6","id17833039-0d7e-47f0-90e6-8cabae94c124","id9a56d8a8-dadd-442e-bc09-732cd619aef7","id178b28cd-99c3-48dc-8d09-1ef71b4cee80","id9a874b64-d0d6-416e-97bc-e9071ed0b16b","id17c1d42c-cb84-4655-a4cd-b54bae17ecaf","id9aa02d46-294b-4e99-b9ee-56edcbfb5275","id17dffffc-65d6-4209-9075-18a441001f0f","id9aa36ac2-8418-4109-a3c1-63ca8742baab","id18c08546-bed3-4ed9-b223-479afe633c8e","id9aabcbd7-a26a-49b1-abf1-02c0dc3adcf2","id1983fade-44f1-41e0-8983-f2ce8d649655","id9af6ed4e-8cdc-4f49-84e9-ba1053b5b3ca","id1987b453-f97d-45c7-9c89-b9a33313d645","id9b132e4f-7e35-4cc5-8711-43ad62b906d0","id19f0cb8c-2e57-4310-967f-a9890f1605db","id9c274536-3ca1-4f0e-93a8-1688074d862f","id1a726069-0b44-42ad-9a20-0ee90b31796b","id9dff11cd-0784-4a60-9c77-43f8957d3dfa","id1adb05d4-e1ca-4ede-81c7-ff24a243abae","id9ecdde9c-2179-4830-bd10-c4fa40f26066","id1b6c184a-5868-4a51-8a82-aa16a7e65126","id9f6be944-83de-42ab-8738-f0022f475e61","id1b703058-e596-45bc-80fe-8b98d545c2e2","id9fb1ba57-2007-4477-b000-2d36f163efd2","id1cc53ade-73c0-492b-a7dc-328269fd0e66","id9fefbe7c-f66a-4940-843e-285cb7b392c1","id1e1cdf72-820f-48a4-8fb6-89575d07cde9","ida06c09a8-0951-49f4-a765-5147bbcae109","id1e308b12-0590-4dae-94d0-a539fcf25df7","ida2ac9937-f351-4d78-9261-264bf6c21e0c","id1ee543d5-b8c0-4f79-8373-6bb6319f2ee2","ida3be445a-a592-4073-815b-47ae41f9c4a0","id1f601832-eee3-48fb-acf5-80c4a454f26e","ida40beabe-496e-4324-a1d1-5a7d8b002aab","id1f6b2aca-7357-40d1-ba7a-99227d9900a2","ida435fdf4-a89a-499e-9006-cbb5693eab26","id2030e10b-a4c2-460f-9e29-90c8a76ea05e","ida482fc4c-ff1d-46a0-a968-7df06422cc4b","id21f6bf91-dfc6-4a59-8b76-88a0f95c7b47","ida5030259-cf9c-4a58-8710-b9da8ee59320","id24118b8f-59fe-482a-aed3-7b1bde1571b7","ida5d3f3fb-6541-419d-b47e-720d438f1bff","id24f21425-b001-4986-aedf-5b4dd851c6ad","ida658becf-3934-49f2-9e74-291dca61d4bb","id2501ee46-8d38-448b-8765-e9c9706cbbe8","ida9d8f797-8a53-4be5-b957-0f1f6d080390","id27988b3e-cde1-4f4b-82a1-6d8ad08db8e9","idad6c9d09-2c03-4786-a72d-dd2aa5f603d4","id2939c03a-6f3f-4f7b-b246-34361baadeb9","idadc853b1-b8bf-488c-a1e2-e95603459b55","id29949b65-2913-4c86-bbd7-b85c5df6f15e","idaecf85cc-058c-46f3-9cdc-3573ac3b8438","id2a84997d-ccee-4f46-bea2-752534f26416","idaf1111af-9384-4e37-9b16-227d38007ebd","id2a971e32-b0e7-4a3f-8c40-62d69c0307d8","idafb00cd4-e5da-4251-86fb-532c3f65324c","id2b14123b-8fcd-402c-9399-4e7c47f20252","idb045f24b-f822-4df9-9ffa-47308edcec8c","id2b76cc9c-c379-43e4-9c3f-56197a3353be","idb094e8b2-8ece-4b36-8025-18073a8b873c","id2cd8ea73-cb3f-4b17-9b97-fa12eb03b85b","idb0cb6ed7-3d9f-4d5e-8df5-edcc921c43b4","id2d00f09f-0074-46c8-b357-b5f33567762a","idb2aac45b-2073-4c7a-adb9-769a4fdcc111","id2d291fa3-f0fd-4bc1-91a6-783863714190","idb46d88e9-8477-4900-8da9-5e7dbded0653","id2d29a4ac-98e7-4663-9dd6-5681bc32ac2e","idb5b99291-507e-4b68-a039-9a0f571f55df","id305f72a8-069b-410b-bbe4-4ecb761c748d","idb5e2cbda-bbfa-4ef8-a9c4-cb978bef9b23","id309005a2-93a8-4566-b8d3-6b9310144266","idb7d6621c-8bce-4a01-8588-9d7761d9372b","id30a82e9f-2a2d-4a66-bdb1-26b881c21d01","idb85b14d7-3b5a-4800-af12-622ec03b9fe5","id30b8f4cd-9245-4496-a8a8-c3e59093bc0a","idb865dec4-f051-4fbe-9405-f832ff2010d7","id31e6e6c8-9197-4927-bfe4-8f6836f963da","idb97bf89a-7a85-4eef-ae7e-f787aead1f0a","id32afdc6c-e899-4795-9e90-c9a192bbc4c1","idba68f2cf-9271-41fd-9655-1fac7681f588","id335d0351-3740-4847-9ff4-be2f078c70f7","idbad34c92-1e64-4187-a0b4-061105fd13ba","id33cf0893-0411-4758-aa90-602bfedf0850","idbbbce1ba-c739-43ba-b9cf-a4f746491ae3","id359f12f9-5c41-48a4-85bc-fd7e307bf7d8","idbce25281-502e-4599-9679-32dc8462ffb1","id3666bc65-8e40-409e-9a1f-41583dd6d978","idbd15f523-45c4-4d7c-a3cc-4fb56abb0e54","id368e23f0-e573-4547-bf5a-14080baf737b","idbd75d8ee-916b-4abc-bf33-bce6a4217076","id374f614a-bcd1-4283-ac6f-fc27ddc79987","idbee9d6c8-948f-4b76-97ca-b985064249ef","id3768d34f-1527-40db-b116-cd80cee5ec3a","idbf339349-062f-4ea9-a0b2-d87d3a21099e","id376dfd27-68e8-4a1a-9c4f-5064279b2a9e","idc0bf9278-9cbb-4361-b4d8-ab172b67e276","id395babb3-3f5d-4e71-a675-af4443f23028","idc2399f5d-dc29-42e9-81f1-1d164d4cc55d","id397d3f69-1453-4057-b177-8723eec923d1","idc2598334-f866-4cc3-93ec-2d2e8b85d319","id39dd61a8-cd17-4df3-ae17-61e239eb00bd","idc2651206-36c5-4f3c-a0d6-3ffc07b27a93","id39eb1db0-2722-4af4-8254-d0e10256c64d","idc3cba26e-38c2-4569-afd6-6b045e604808","id3acf7438-fd65-442c-8a30-68b6714537f3","idc435627c-159d-4a6d-a819-30abac24bf4d","id3c037acf-f453-4513-a6dd-129163ddde2a","idc4d1e105-28b7-48df-abef-fbe09782fdb2","id3d1f4059-2220-45b4-a4d2-b14f76cec96a","idc7df3466-b9a7-4818-883b-d0cd08483570","id3d2aa654-1b5f-4eb4-a1c2-af31f5760069","idc88e3901-f3ee-436a-ba80-143c8eab7d69","id3e5f451a-5882-4914-ae1a-95c898c2bcca","idc8febeef-8e7e-459e-88cb-8086692dc559","id3ec3bcb3-492c-4942-be02-76a67f57e38f","idc9e28934-8379-4511-817c-d787f2c4ca3a","id3ed614e7-f356-4d87-985b-d3bbbae3bb40","idca554128-da9d-4f37-9560-ca083509e01d","id3fc3755d-a3f8-4e2c-813f-ff124f2a75c1","idcb99fb14-641c-4fb3-bf1d-17efb8cd982a","id419deaac-ea45-4bdd-9fa0-b5cd8429b44f","idcce62116-3dcd-400d-ada0-9839ca02466e","id42251e02-b687-4bf1-a1b6-e3dd978542e4","idcede5274-279d-4021-8aca-5ecc4bf94d66","id422a46b2-a67c-4a7e-923f-9b651ced96f8","idd3164236-c14a-4230-b527-ecd1a3992f02","id44bec761-b603-49c0-8634-f6bfe0319bb1","idd32a4219-4cca-470f-b6b3-11447d9f5e0b","id45bdcfd6-1e3f-4be8-b843-ae949e8e43eb","idd39df936-4bc6-404d-9e1d-f6416cc5a34c","id45d58088-d3f4-4eb0-b75c-bcb98c0981dc","idd3d545b3-457f-4389-821f-704cb24aff7f","id49717f75-0f2d-4e1c-9a12-f1cd7877b80a","idd4bc755a-2585-4529-ae36-7e1d88bdecfe","id4b930a10-4b12-4428-84f9-3255b4a3bc4f","idd52c3cb6-921e-4f26-b49c-e6e4f11f9c06","id4c18d9cf-4af4-4a86-8b1c-f78795fbbd7e","idd58b535a-6c95-43bb-ad88-0600e6447537","id4c42dc4e-66b7-40bf-9fbe-b92543248198","idd680df09-368e-42b5-b540-45c41ed31042","id4ca13c92-84b4-4edf-842a-b20b7e713415","idd6f7afc0-1558-43ad-acb1-2b5311ed2264","id4d02010e-9e69-40d3-9bf0-3d4510aa0614","idd77e4dbc-b239-4742-8cb7-efd427010d13","id4da999a0-ef41-4a0b-b1d1-446b39cc855a","idd89b1fd6-bef4-4803-8ed3-3b442be600b6","id4e6edfe6-adcb-4c12-8ff4-38a79f5887e8","idd9d7f0d4-fa64-48f2-8037-52cbf0b34ad8","id4e9d13da-438d-4b5e-8e86-67432b5e9471","idd9dc3b59-613d-469e-8b4f-6c5a557eb26a","id4eac2c98-86d2-4ee6-a1d3-157d013c78dc","idda93a143-5799-4856-a5f4-3ff3b9284311","id4fc8e011-4433-4537-b03f-457a3a70240f","iddb1e3104-8bfd-4353-854a-c478ba1f7d07","id5017c115-39b8-4c5c-8be6-2cb298d7c641","iddc6ceece-efaa-4b9d-9e12-d35863d7a902","id501c987e-d1eb-48a9-89eb-72a5062c90b4","iddca004bf-9c14-45a1-b186-ab3759f6c7fe","id52292ffc-0902-4d97-b461-20723987a177","idde287a25-981d-49f8-81af-fe0300e4cd31","id53886143-c1c6-40e9-88e6-e4e5e0271fc8","idde7b7cac-f094-4d59-8651-e991e34ea093","id53d9a130-4a3a-4491-85cd-9a44415a1632","iddf576520-a6b6-4c9b-8d06-3f59cc5342fd","id556fcbc8-172a-4af1-8822-ae036e8d68e8","iddf5ab6cc-6f68-4b6b-95e2-954c6b57ba9c","id5613f5f6-086d-470b-8f77-6dcb7f8625b7","iddfd30f94-19e3-4838-95df-14a3763419af","id56a82a56-0241-4d3d-9de2-696b0c36df91","ide17a6048-7a72-42c0-ad3f-97cbff02bc9f","id577847b9-f9b6-4954-868f-3d5ab2b4f694","ide21964b8-47c7-4a5f-bfa7-7206222883d2","id58d34254-4f5b-40a4-9e9f-7160062fb2a4","ide303b8c4-6428-4134-bb9d-a95dc1d8e35e","id5c127332-5ca0-45f1-a5ac-4876ad94e491","ide3711a9b-6d4c-44df-bbab-0a675046a5df","id5c22c8d6-5805-455c-8e1c-e6b8006738e4","ide43d3769-e25f-4fb7-9080-ea26defaf094","id5c984433-33cf-42fc-b3ba-511efcdcab19","ide4b08e66-0772-4956-8670-e5c6553eecbe","id5d2f1a0f-a01b-48c2-b749-1ca710aa2bb4","ide56491af-718f-4ca1-b3cd-0e5b3154697b","id5d54c742-5a8e-4c40-8d62-95e75e210ab8","ide5e4aa9f-6015-4d62-a088-69920d7f2500","id5d8ed961-a012-4b6b-bffb-9468d3953646","ide6827400-0d95-46b0-8874-6ce9e9d5011b","id5dd423e8-feaa-4568-a750-500948c41d6c","ide6b72c24-1607-43b9-8b8a-7bf83eea5895","id5e18b17d-4626-4b6d-8ac6-e560cee0376c","ide733f289-a0ec-478c-aad9-1378ad26ddfa","id5ed024e8-d05e-4c65-9441-eda9930ccc82","ide7d1b042-5b1c-4c29-b501-a5ef4810ecd2","id5fce954c-11b9-41f4-9471-4b5bdc0efa7f","ide7d49212-e9cb-4cd7-923a-bdfc82a4a150","id5fd9552a-c742-4388-940d-295d1107ae00","ide8ac39e6-377a-412c-b18d-f6c4fc654ec6","id5ffa2ec6-2d94-4b09-8fb9-3591cf38fa4f","ide9483296-cb91-497a-b955-39a3c3289dac","id60df7543-6da5-4c75-943b-5800c1e08234","ide95e9e7d-a1ed-46eb-9cf8-1d355bf0577f","id6101ffe6-2ef6-4256-9d6b-a7c545836995","idea3b5da2-6a12-400c-bf0f-e442f5ec1132","id6176621d-8534-4fba-a0e7-74a5a78245a9","ideadf8482-e60e-4307-adb7-d5c3b9fa6cae","id61eb08d0-ea01-45e7-884b-ac3d6070172a","ideaf674dd-9909-4879-b9a3-74253169eddd","id62379be5-13f0-474b-94d3-6f944ec4ee96","ideb59cf69-1997-41b9-bf69-69ad7da292a1","id648252da-dd86-48d1-ae75-4257e3142b0b","idebcba7f2-ce13-4bae-97cd-91a6b1dcd465","id64bb5550-2735-4401-a0db-58ec1020a32d","idee850103-6314-41a3-9734-8b520de95b40","id659294b9-ded9-498a-bc00-c1d4d456aa4a","ideeb9d147-608d-4692-8adf-2f601d23a8ff","id65cac997-4d39-4501-85ec-4fcb328a8eb5","idef4cbd38-bc79-4d60-a715-647edd2ebe9e","id65ddccd2-70a7-48b3-88e1-9d6fd1dcd11b","idf030e410-2585-461e-ad23-9fdf026ac06c","id66c92b9e-de3c-4d1a-bb69-46ffbc6caf33","idf05d314c-5ec5-4e2a-b785-9a702716f111","id66e305ca-8a16-44f0-b837-85e529b7ad55","idf063ddbb-1668-40df-9ec2-b0a23ca2c389","id670c1cbc-4494-45f1-bb8b-18db82d4f7e0","idf0787165-6f58-4d67-b510-928eea2c4882","id6909302d-358a-46cd-8052-eedb09c83c83","idf0c353fd-947c-41e2-b643-3ecc0d69796c","id69d56f2d-6924-409b-9d1e-c8d69b400270","idf33aa93f-b32d-4c4f-9a71-1e5ab69c6899","id6b74c07a-9932-4a80-8296-9db08d0c6964","idf3ba1cf8-ac67-4da8-96f3-8680dc55f247","id6c28f086-6a25-40b6-93eb-bba0014acda6","idf435ec00-0db8-46b4-96bf-890a9d931df3","id6d10d4ee-6331-4bba-93bc-a7b64cc0b22a","idf45210d3-9e66-4f5e-bef1-5ee5547cc893","id6e6962d8-34b8-431c-8220-42b0b92a410b","idf5033d52-ec36-4e67-88d8-b6d898b81b2f","id6fa2a667-9c36-4526-8a58-1975e863a806","idf51af6db-2655-47fb-9ffc-83a503a728ea","id700e91bb-d675-41b2-bbbd-935767c7b447","idf6e5979d-d6a9-4408-bd1b-1b37d1064e83","id70741b73-9683-42bd-87fc-6a021a70103b","idf7e06a10-74ca-4f52-a944-b3f22c016998","id70939245-41e0-4845-a473-5fef719b9828","idf86fa219-34e9-4812-8edf-2b66d7219779","id74139255-a635-4c87-814d-3dd04ed630a8","idf8c0a69c-df53-4a0a-a30e-4424c0ce3dd5","id748e38b1-2ead-4a0a-8881-d640617b856b","idf9ceffc0-d544-418d-b4a9-bd3c84e37026","id74a0264d-1d31-430d-9a88-e7334c8aa96c","idf9cf605f-f287-4eba-8e8f-8ba47e14ecc4","id74f31744-0aed-4633-beec-4e203315f0b7","idfa1f1cfc-9da0-488b-93fe-60aaff99c171","id750770db-53cf-4021-a8c2-5e0ddb97fef5","idfafd4576-0709-49e1-926f-ad1ee09907ec","id75113445-d2d6-44a0-866c-c9175e6d214b","idfafd6f5b-1d76-4537-bd1c-e0bd7b4e2166","id7664c241-b369-4d76-9cd1-7f2bdc1d267f","idfb66b7f5-22ca-4d4d-97c2-30e6edfff781","id78048432-2de1-498c-9cf1-9c5bd43daa39","idfb9bafa5-7133-4955-8156-4eb6763dc8e1","id78412f82-cc55-43d5-96ca-1d02cf957725","idfc3b7596-f515-446f-81db-fed0154ca2c5","id785803bb-42a8-4afb-b7b4-58c114894408","idfef9c64f-5959-4da0-aaa2-66b56fc7b4c3","id790ff5db-b7f3-4946-aaed-305f66b1dd6a","id7941d1c7-5319-4e8d-a15f-60b61cb68a9b"]
#errored "id01eef340-598c-4205-a990-cec190ac2ca5"
if __name__ == '__main__':
    for pid in pids:
        runWithUnzip(WGSDir, outputDir, pid)
#         id = filepath.split("/")[0]
#         print id
#         outputfile="varscan.somatic.snp.coverage.summary"
#         try: run(os.path.join(WGSDir, filepath), os.path.join(outputDir, outputfile), id)
#         except: print id, "failed"
        
        
#         run("/Users/cmelton/Documents/AptanaStudio3Workspace/CODES/TestData/Somatic_SNVs_Indels/cp2.varscan.snp",
#         "/Users/cmelton/Documents/AptanaStudio3Workspace/CODES/TestData/Somatic_SNVs_Indels/cp2.varscan.coverage.summary")
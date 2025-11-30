import PyPDF2

# 打开PDF文件
with open(r'd:\2026YJ\论文写作\图注.pdf', 'rb') as file:
    # 创建PDF阅读器对象
    reader = PyPDF2.PdfReader(file)
    
    # 获取PDF的页数
    num_pages = len(reader.pages)
    print(f"PDF总页数: {num_pages}\n")
    
    # 逐页读取内容
    for page_num in range(num_pages):
        page = reader.pages[page_num]
        text = page.extract_text()
        print(f"第 {page_num + 1} 页内容:")
        print("=" * 50)
        print(text)
        print("\n" + "=" * 50 + "\n")
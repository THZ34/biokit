# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity


def text_similarity(input_list, ref_list):
    """匹配两个列表中相似度最高的文本

    :param input_list:
    :param ref_list:
    :return:
    """
    vectorizer = TfidfVectorizer()

    all_texts = input_list + ref_list
    tfidf_matrix = vectorizer.fit_transform(all_texts)

    # 分离a和b的TF-IDF矩阵
    a_tfidf = tfidf_matrix[:len(input_list)]
    b_tfidf = tfidf_matrix[len(input_list):]

    # 用于存储最终结果的字典
    result = {}

    # 对于a列表中的每个文本，找到b列表中最相似的文本
    for i, text_a in enumerate(input_list):
        # 计算a[i]与b中每个文本的余弦相似度
        similarities = cosine_similarity(a_tfidf[i], b_tfidf)

        # 找到b中与a[i]相似度最高的文本及其相似度
        best_match_index = similarities.argmax()
        best_match = ref_list[best_match_index]
        best_similarity = similarities[0][best_match_index]

        # 将结果添加到字典中
        result[text_a] = (best_match, best_similarity)

    # 输出结果
    return result

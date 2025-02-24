# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity


def text_similarity(input_list, ref_list):
    """ƥ�������б������ƶ���ߵ��ı�

    :param input_list:
    :param ref_list:
    :return:
    """
    vectorizer = TfidfVectorizer()

    all_texts = input_list + ref_list
    tfidf_matrix = vectorizer.fit_transform(all_texts)

    # ����a��b��TF-IDF����
    a_tfidf = tfidf_matrix[:len(input_list)]
    b_tfidf = tfidf_matrix[len(input_list):]

    # ���ڴ洢���ս�����ֵ�
    result = {}

    # ����a�б��е�ÿ���ı����ҵ�b�б��������Ƶ��ı�
    for i, text_a in enumerate(input_list):
        # ����a[i]��b��ÿ���ı����������ƶ�
        similarities = cosine_similarity(a_tfidf[i], b_tfidf)

        # �ҵ�b����a[i]���ƶ���ߵ��ı��������ƶ�
        best_match_index = similarities.argmax()
        best_match = ref_list[best_match_index]
        best_similarity = similarities[0][best_match_index]

        # �������ӵ��ֵ���
        result[text_a] = (best_match, best_similarity)

    # ������
    return result

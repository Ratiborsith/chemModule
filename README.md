# Программный модуль анализа химических формул в патентных документах

## Достигнутая цель
Разработан программный модуль, способный распознавать и классифицировать информацию, а затем преобразовывать полученные данные в единый формат Fingerprint для последующего сравнения.
<br>
Область применения данной разработки – исследование химических патентов, анализ, распознавание и классификация химических представлений в них. 

## Решенные задачи
Для достижения данной цели были осуществлены следующие задачи задачи:
1. Проведен анализ предметной области, рассмотреть системы-аналоги, методы решения задачи распознавания графических и текстовых представлений химических формул, преобразования различных типов формул (IUPAC, SMILES, Inchi, MOL, структурная формула).
2. Разработан алгоритм распознавания графических представлений химических формул.
3. Разработан алгоритм распознавания текстовых представлений химических формул.
4. Разработан алгоритм преобразований различных представлений химических формул.
5. Спроектирован программный модуль.
6. Программно реализован модуль и проверена его эффективность.

## Установка
Чтобы запустить этот программный модуль анализа химических формул в патентных документах, убедитесь, что на вашей машине(MacOs 14 Sonoma или Windows 10+) установлен Python 3.10 или новее. Клонируйте этот репозиторий, а затем установите необходимые зависимости:
```bash
git clone https://github.com/Ratiborsith/chemModule
```

## Использование
Для того, чтобы запустить графический интерфейс программного модуля введите следующую команду в терминале:
```bash
python app.py
```
После запуска вас встретит простой и интуитивно понятный интерфейс. Вы сможете осуществлять следующие функции:
1.	 Сравнение по коэффициенту Tanimoto химических формул разных представлений друг с другом. Доступны для сравнения следующие химические представления: InChi, MOLFILE, SMILES, структурная формула, IUPAC
2.	Сравнение патентных документов. Сравниваются по парам все химические формулы, которые присутствуют в патенте.
3.	Парсинг структурных и IUPAC представлений из Google Patents
4.	Визуализация представленной информации в удобном для пользователя формате. 


## Технологии
#### • 	СУБД: SQLite3 для хранения данных.
#### •	Библиотеки для химического анализа: RDKit, Cirpy.
#### •	Веб-сервер: Flask.
#### •	Средства аутентификации и авторизации: Flask-Login.
#### •	Скрипты интерфейса: JavaScript, ajax
#### •	Внешний вид веб-интерфейса: html, css, scss. Для работы с формами: flask wtforms.
#### •	Библиотеки для подмодуля парсинга и краулинга: Selenium, BeautifulSoup.
#### •	Используемая модель для конвертации структурных формул: OSRA.

## Входные данные
#### •	Ссылки на химические патенты; 
#### •	Информация о патентах; 
#### •	Данные о требуемых пользователем форматов соединений;
#### •	Для анализа IUPAC, MOLFILE и структурных представлений занесенные соответствующие представления(или пути к ним) в базу данных.

## Выходные данные
#### •	Результат сходства двух химических формул на странице сравнения химических формул: 
По запросу пользователя модуль выводит результат сравнения двух формул разных представлений. Выводится коэффициент Tanimoto и структурные(графические) представления двух представлений.
#### •	Список формул, схожих в выбраных патентах:
По запросу пользователя модуль выводит результат сравнения двух патентов. Выводятся в порядке убывания по коэффициенту Tanimoto все формулы в данных патентах, сравненные друг  с другом.
#### •	Информация о патентах:
Модуль может предоставлять дополнительную информацию о патентах, включая название, авторов, дату публикации и другие метаданные.
#### •	Выходные данные в единый формат SMILES на странице сравнения патентов:
После анализа и сравнения химических формул, модуль преобразует полученные данные в единый формат SMILES для удобства дальнейшего использования и сравнения.

## Архитектура
<p align="center">
<img width="600" alt="image" src="https://github.com/Ratiborsith/chemModule/assets/85187788/2aceb93e-ce1f-4e14-b350-1a5c366853e4"> <br>
</p>
1. – обращение к базе Google Patents по ссылке патента;
2. – получение из Google Patents кода страницы патента;
3. – вызов блока парсинга формул;
4. – передача в БД спарсенных химических представлений IUPAC и структурной формулы;
5. – получение из БД списка патентов;
6. – получение из БД списка патентов и содержащихся в каждом выбранном патенте химических соединений;
7. – вызов блока создания молекул RDKit;
8. – вызов блока получения fingerprint;
9. – Обращение к подмодулю сравнения. Вызов блока сравнения формул по коэффициенту Tanimoto для двух fingerprint;
10. – вызов блока сравнения патентов;
11. – обращение к подмодулю визуализации. вызов блока результатов рассчета коэффициента Tanimoto;
12. – обращение к БД для проверки авторизации;
13. – получение из БД пароля для указанного по логину пользователя. Сравнение через Flask-Login с указанными;
14. – получение из БД списка патентов.


## Как выглядит программный модуль?
<p align="center">
<img width="600" alt="image" src="https://github.com/Ratiborsith/chemModule/assets/85187788/f905b260-0e09-429a-aa2c-77b8abb8a53c">

<br>
Рисунок 1 – Авторизация
<br>
<img width="600" alt="image" src="https://github.com/Ratiborsith/chemModule/assets/85187788/90851b3c-e77c-41c0-ae12-995c325b718b">
<br>
Рисунок 2 – Список патентов
<br>
<img width="600" alt="image" src="https://github.com/Ratiborsith/chemModule/assets/85187788/ebff24fb-9354-4472-8e8a-7eac4cfc9bbe">
<br>Рисунок 3 – Парсер<br>
<img width="600" alt="image" src="https://github.com/Ratiborsith/chemModule/assets/85187788/14a3efa6-82b0-4bac-94f3-bfc8662c5863">
<br>Рисунок 4 – Сравнение двух химических формул<br>
<img width="600" alt="image" src="https://github.com/Ratiborsith/chemModule/assets/85187788/c0261af6-4552-46d4-bda2-25f27636862c">
<br>Рисунок 5 – Сравнение двух патентов<br>
</p>





